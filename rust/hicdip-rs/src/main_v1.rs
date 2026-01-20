use clap::Parser;
use std::cmp::min;
use hdf5::{file, types, Error, File, Group, Result,Dataset};
use std::collections::HashMap;
use polars::prelude::*;
use polars::io::csv::CsvWriter;
use hdf5::types::{VarLenArray, FixedAscii, TypeDescriptor};
use ndarray::Array2;
use rand::seq::index::sample;
use rand::thread_rng;
use std::fs::File as STDFILE;
use std::path::{Path, PathBuf};

//evcxr

//CLI argument parser
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args{
    #[arg(long)]
    hap1: String,
    #[arg(long)]
    hap2: String,
    #[arg(long)]
    reference: String,
    #[arg(long)]
    chr_sizes: String,
    #[arg(long, default_value = "100000")]
    bin_size: u32,
    #[arg(long)]
    chr_x: String,
    #[arg(long)]
    chr_y: String,
    #[arg(long)]
    skip_diags: u16,
    #[arg(long, default_value = "1000")]
    n_permutations: u32,
}

pub fn get_n_chroms(file: &File) -> hdf5::Result<usize> {
    let group = file.group("chroms")?;
    let dataset = group.dataset("name")?;
    Ok(dataset.size())
}

pub fn read_cool(path: &String) -> hdf5::Result<File> { 
    let cool = File::open(path)?;
    Ok(cool)
}

pub fn get_bins_group(file : &File, resolution:&u32) -> hdf5::Result<Group>{
    let bins = file.group(&format!("resolutions/{}/bins",&resolution))
    .expect("No bins group present");
    Ok(bins)
}

pub fn get_chr_start_end(group : &Group) -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<u32>)>{
    let chrom  = group.dataset("chrom")?.read_1d::<u32>()?.to_vec(); // Get the bin1_id dataset from hap1
    let start  = group.dataset("start")?.read_1d::<u32>()?.to_vec(); // Get the bin2_id dataset from hap1
    let end = group.dataset("end")?.read_1d::<u32>()?.to_vec(); // Get the count dataset from hap1
    Ok((chrom,start,end))
}

pub fn get_pixels(file: &File, resolution: &u32) -> hdf5::Result<(Vec<u32>, Vec<u32>, Vec<f32>)> {
    let pixels_group = file.group(&format!("resolutions/{}/pixels", resolution))?;

    let bin1_id = pixels_group.dataset("bin1_id")?.read_1d::<u32>()?.to_vec();
    let bin2_id = pixels_group.dataset("bin2_id")?.read_1d::<u32>()?.to_vec();
    let counts = pixels_group.dataset("count")?.read_1d::<f32>()?.to_vec(); // could also be u32 or f64

    Ok((bin1_id, bin2_id, counts))
}

//Vec<u32>,Vec<u32>,Vec<u32>
pub fn kr_norm(file:&File,resolution:&u32) -> hdf5::Result<(Vec<u32>,Vec<u32>,Vec<f32>)> {
    let (bin1, bin2, raw_count) = get_pixels(file, resolution)?;
    let bins_group = get_bins_group(file, resolution)?; 
    let kr_w = bins_group.dataset("KR")?.read_1d::<u32>()?.to_vec(); // gives a vector with the weights 
    let raw_iterator = raw_count.iter().zip(bin1.iter().zip(bin2.iter())); // (&raw,(&bin1,&bin2))
    let norm_counts: Vec<f32> = raw_iterator.map(|(&count,(&i,&j))|{
        // now we are inside the raw counts for the matrix and because they are the same dimension we can use i,j to retrieve3 the coefficients from KR
        let w_i = kr_w.get(i as usize).copied().unwrap_or(0.0 as u32);
        let w_j= kr_w.get(j as usize).copied().unwrap_or(0.0 as u32);
        if w_i > 0 && w_j > 0 {
            count / (w_i as f32 * w_j as f32)
        } else {
            0.0
        }
    })
    .collect();

    Ok((bin1,bin2,norm_counts))

}

pub fn get_chromosomes(file: &File,resolution:&u32) -> hdf5::Result<Vec<String>> {
    let chr_names = file.group(&format!("/resolutions/{}/chroms", &resolution))?;
    let dataset_chr = chr_names.dataset("name")?;
    let raw_names = dataset_chr.read_1d::<FixedAscii<32>>()?;
    let chrom_names: Vec<String> = raw_names
        .iter()
        .map(|s| s.as_str().trim_end_matches('\0').to_string())
        .collect();
    Ok(chrom_names)
}



pub fn check_df_index(target1: u32,
                    target2: u32,
                    dataframe: DataFrame,) -> Result<(),Error>{

    let filtered = dataframe
                .lazy().filter(
                        col("bin1").eq(lit(target1))
                            .and(col("bin2").eq(lit(target2))),
                        ).collect().map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;
    println!("{:?}",filtered);
    Ok(())
}



pub fn get_hic_matrix(
    file: &File,
    resolution: &u32,
    chr_group: &Group,
    chr_x: &str,
    chr_y: &str,
) -> hdf5::Result<DataFrame> {
    let chr_names= get_chromosomes(file, resolution)?; // names of chromosomes
    ////////////
    /// // whats the index of the chromosome 
    let chr_x_index = chr_names.iter()
                            .position(|c| c == chr_x || c == &format!("chr{}",chr_x))
                            .unwrap(); //0 based     // Get the index of the chromosome to read
    let chr_y_index = chr_names.iter()
                            .position(|c| c == chr_y || c == &format!("chr{}",chr_y))
                            .unwrap();

    // Just get all the bins 
    let bin_ids = chr_group.dataset("chrom")?.read_1d::<u32>()?.to_vec();
    // enumerate wraps the iterator composed of (i,&c) wherre i is the  index of the bin and c is the chromosome index
    let x_indeces: Vec<usize> = bin_ids.iter().enumerate()
                        .filter_map(|(i,&c)| if c == chr_x_index as u32 {Some(i)} else {None})
                        .collect(); 
    let y_indeces: Vec<usize> = bin_ids.iter().enumerate()
                    .filter_map(|(i,&c)| if c == chr_y_index as u32 {Some(i)} else {None})
                    .collect(); 

    // then we get the indexes that correspond to our chromosome 
    let (min_x,max_x) = (*x_indeces.first().unwrap(),*x_indeces.last().unwrap());
    let (min_y,max_y) = (*y_indeces.first().unwrap(),*y_indeces.last().unwrap());

    // Read the hic matrix pixels 
    let (bin1, bin2, count) = get_pixels(file, resolution)?;
       
    let mut contact_map = HashMap::new();
    for i in 0..bin1.len(){
       let b1 = bin1[i] as usize;
       let b2 = bin2[i] as usize;
       let val = count[i];

       if (min_x..=max_x).contains(&b1) && (min_y..=max_y).contains(&b2){
           contact_map.insert(((b1-min_x) as u32, (b2-min_y) as u32),val);
       } else if (min_y..=max_y).contains(&b1) && (min_x..=max_x).contains(&b2){
           contact_map.insert(((b1-min_y) as u32, (b2-min_x) as u32),val);
       }
     }

     let x_len = max_x - min_x +1;
     let y_len = max_y - min_y + 1;

     let mut bin1_vec = Vec::with_capacity(x_len*y_len);
     let mut bin2_vec = Vec::with_capacity(x_len*y_len);
     let mut count_vec = Vec::with_capacity(x_len*y_len);

     for i in 0..(x_len as u32) {
       for j in 0..(y_len as u32) {
           bin1_vec.push(i);
           bin2_vec.push(j);
           count_vec.push(*contact_map.get(&(i,j)).unwrap_or(&0.0))
       }
     }

     let df = DataFrame::new(vec![
       Series::new("bin1",bin1_vec),
       Series::new("bin2",bin2_vec),
       Series::new("count",count_vec),
     ]).map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;

    Ok(df)
}


pub fn get_hic_matrix_kr(
    file: &File,
    resolution: &u32,
    chr_group: &Group,
    chr_x: &str,
    chr_y: &str,
) -> hdf5::Result<Array2<f32>> {
    let chr_names = get_chromosomes(file, resolution)?;
    
    let kr_w = file
        .group(&format!("resolutions/{}/", resolution))?
        .dataset("bins/KR")?
        .read_1d::<f32>()?
        .to_vec();

    let chr_x_index = chr_names
        .iter()
        .position(|c| c == chr_x || c == &format!("chr{}", chr_x))
        .unwrap();
    let chr_y_index = chr_names
        .iter()
        .position(|c| c == chr_y || c == &format!("chr{}", chr_y))
        .unwrap();

    let bin_ids = chr_group.dataset("chrom")?.read_1d::<u32>()?.to_vec();

    let x_indices: Vec<usize> = bin_ids
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| if c == chr_x_index as u32 { Some(i) } else { None })
        .collect();
    let y_indices: Vec<usize> = bin_ids
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| if c == chr_y_index as u32 { Some(i) } else { None })
        .collect();

    let (min_x, max_x) = (*x_indices.first().unwrap(), *x_indices.last().unwrap());
    let (min_y, max_y) = (*y_indices.first().unwrap(), *y_indices.last().unwrap());

    let x_len = max_x - min_x + 1;
    let y_len = max_y - min_y + 1;

    // Initialize dense matrix
    let mut matrix = Array2::<f32>::zeros((x_len, y_len));

    let (bin1, bin2, count) = get_pixels(file, resolution)?;

    for i in 0..bin1.len() {
        let b1 = bin1[i] as usize;
        let b2 = bin2[i] as usize;
        let val = count[i];

        let w_1 = kr_w.get(b1).copied().unwrap_or(0.0);
        let w_2 = kr_w.get(b2).copied().unwrap_or(0.0);

        if w_1 > 0.0 && w_2 > 0.0 {
            let normalized_val = val / (w_1 * w_2);
            if (min_x..=max_x).contains(&b1) && (min_y..=max_y).contains(&b2) {
                matrix[(b1 - min_x, b2 - min_y)] = normalized_val;
            } else if (min_y..=max_y).contains(&b1) && (min_x..=max_x).contains(&b2) {
                matrix[(b2 - min_x, b1 - min_y)] = normalized_val;
            }
        }
    }

    Ok(matrix)
}

pub fn get_hic_df_kr(
    file: &File,
    resolution: &u32,
    chr_group: &Group,
    chr_x: &str,
    chr_y: &str,
) -> hdf5::Result<DataFrame> {
    let chr_names= get_chromosomes(file, resolution)?; // names of chromosomes
    // for the whole matrix
    let kr_w = file.group(&format!("resolutions/{}/", resolution))?
               .dataset("bins/KR")?
               .read_1d::<f32>()?
               .to_vec();
    ////////////
    /// // whats the index of the chromosome 
    let chr_x_index = chr_names.iter()
        .position(|c| c == chr_x || c == &format!("chr{}",chr_x))
        .unwrap(); //0 based     // Get the index of the chromosome to read
    let chr_y_index = chr_names.iter()
        .position(|c| c == chr_y || c == &format!("chr{}",chr_y))
        .unwrap();
    // Just get all the bins 
    let bin_ids = chr_group.dataset("chrom")?.read_1d::<u32>()?.to_vec(); // the value at each position is the index of the chromosome that bin belongs to 
    // let kr_w = chr_group.dataset("KR")?.read_1d::<u32>()?.to_vec();
    // enumerate wraps the iterator composed of (i,&c) wherre i is the  index of the bin and c is the chromosome index
    let x_indeces: Vec<usize> = bin_ids.iter().enumerate()
                        .filter_map(|(i,&c)| if c == chr_x_index as u32 {Some(i)} else {None}) // because we have the index and the values of the vector are the indexes we can compare
                        .collect(); 
    let y_indeces: Vec<usize> = bin_ids.iter().enumerate()
                    .filter_map(|(i,&c)| if c == chr_y_index as u32 {Some(i)} else {None})
                    .collect(); // global bin numbers 
    
    // Kr for chromosome
    let kr_w_x: Vec<f32>  = x_indeces.iter().map(|&i| kr_w[i] as f32).collect(); // with i coming from any of our indexes retrieve the value from kr 
    let kr_w_y: Vec<f32>  = y_indeces.iter().map(|&i| kr_w[i] as f32).collect();

    // then we get the indexes that correspond to our chromosome 
    let (min_x,max_x) = (*x_indeces.first().unwrap(),*x_indeces.last().unwrap());
    let (min_y,max_y) = (*y_indeces.first().unwrap(),*y_indeces.last().unwrap());

    // Read the hic matrix pixels, all hic 
    let (bin1, bin2, count) = get_pixels(file, resolution)?;
    
    let mut contact_map = HashMap::new();
    for i in 0..bin1.len(){
       let b1 = bin1[i] as usize; // get the global bin 
       let b2 = bin2[i] as usize;
       let val = count[i];

       let w_1 = kr_w.get(b1).copied().unwrap_or(0.0);
       let w_2= kr_w.get(b2).copied().unwrap_or(0.0);

       if w_1 > 0.0 && w_2 > 0.0 {
        let normalized_v = val / (w_1 as f32 * w_2 as f32);
        if (min_x..=max_x).contains(&b1) && (min_y..=max_y).contains(&b2){
            contact_map.insert(((b1-min_x) as u32, (b2-min_y) as u32),normalized_v);
        } else if (min_y..=max_y).contains(&b1) && (min_x..=max_x).contains(&b2){
            contact_map.insert(((b1-min_y) as u32, (b2-min_x) as u32),normalized_v);
        }

       }

    
     }
     /// til here we have a matrix
     
     let x_len = max_x - min_x +1;
     let y_len = max_y - min_y + 1;

     let mut bin1_vec = Vec::with_capacity(x_len*y_len);
     let mut bin2_vec = Vec::with_capacity(x_len*y_len);
     let mut count_vec = Vec::with_capacity(x_len*y_len);

     for i in 0..(x_len as u32) {
       for j in 0..(y_len as u32) {
           bin1_vec.push(i);
           bin2_vec.push(j);
           count_vec.push(*contact_map.get(&(i,j)).unwrap_or(&0.0))
       }
     }

     let df = DataFrame::new(vec![
       Series::new("bin1",bin1_vec),
       Series::new("bin2",bin2_vec),
       Series::new("count",count_vec),
     ]).map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;

    Ok(df)
}

pub fn dense2tag(matrix:&Array2<f32>) -> (Vec<(usize, usize)>, usize){
    let mut tag_mat = Vec::new();

    let shape = matrix.shape();
    let rows = shape[0];
    let cols = shape[1];

    for i in 0..rows{
        for j in i..cols { // Upper
            let count = matrix[(i,j)].ceil() as u32;
            for _ in 0..count{
                tag_mat.push((i,j)); // each contact between bins i,j is represented as a row
            }
        }
    }
    let tag_len = tag_mat.len(); // because we generated N numbers of rows per count

    (tag_mat,tag_len) 

}

pub fn tag2dense(tag: &[(usize, usize)], nsize: usize) -> Array2<u32> {
    let mut count_map: HashMap<(usize, usize), u32> = HashMap::new();

    for &(i, j) in tag {// count occurrences of each pair
        *count_map.entry((i, j)).or_insert(0) += 1;
    }


    let mut dense = Array2::<u32>::zeros((nsize, nsize)); //  zero matrix

    // 
    for (&(i, j), &val) in &count_map {
        dense[(i, j)] = val;
        if i != j {
            dense[(j, i)] = val;  // Make symmetric
        }
    }

    dense
}

fn random_sampling(tag_len:usize,n_reads:usize) -> Vec<usize>{
    let mut rng = thread_rng();
    let sample_idx = sample(&mut rng,tag_len,n_reads)
                    .into_vec();
    sample_idx

}

pub fn downsampling(matrix:&Array2<f32>,n_reads:usize)-> Array2<u32>{

    let (tag_mat,tag_len) = dense2tag(&matrix);
    let sample_idx = random_sampling(tag_len,n_reads); // from tag_len select n_reads amount 
    let sample_tag: Vec<(usize,usize)> = sample_idx.iter().map(|&i| tag_mat[i]).collect(); 

    let rows = matrix.shape()[0];

    let down_matrix = tag2dense(&sample_tag,rows);
    down_matrix
    

}

fn array2dataframe(matrix:&Array2<u32>) -> Result<DataFrame> {
    let shape = matrix.shape();
    let rows = shape[0];
    let cols = shape[1];

    let mut bin1_vec = Vec::new();
    let mut bin2_vec = Vec::new();
    let mut count_vec = Vec::new();

    for i in 0..rows {
      for j in 0..cols {
          bin1_vec.push(i as u32);
          bin2_vec.push(j as u32);
          count_vec.push(matrix[(i,j)]);
      }
    }

    let df = DataFrame::new(vec![
      Series::new("bin1",bin1_vec),
      Series::new("bin2",bin2_vec),
      Series::new("count",count_vec),
    ]).map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;

    Ok(df)
    
}


fn main() -> Result<()>{
    let args = Args::parse();
    let hap1 = read_cool(&args.hap1)?;
    let hap2 = read_cool(&args.hap2)?;
    let reference = read_cool(&args.reference)?;
    let n_permutations = args.n_permutations;

    let hap1_group = get_bins_group(&hap1,&args.bin_size)?;
    let hap2_group = get_bins_group(&hap2,&args.bin_size)?;
    let reference_group = get_bins_group(&reference,&args.bin_size)?;

    let contacts_hap1 = get_hic_df_kr(&hap1, &args.bin_size, &hap1_group, &args.chr_x,  &args.chr_y)?;
    let contacts_hap2 = get_hic_df_kr(&hap2, &args.bin_size, &hap2_group, &args.chr_x,  &args.chr_y)?;
    let contacts_reference = get_hic_matrix_kr(&reference, &args.bin_size, &reference_group, &args.chr_x,  &args.chr_y)?;

    
    let mut joined_df = contacts_hap1.join(
        &contacts_hap2,
        ["bin1", "bin2"],
        ["bin1", "bin2"],
        JoinArgs::new(JoinType::Inner),
    )
    .map_err(|e| hdf5::Error::Internal(format!("Polars join error: {}", e)))?;

    let count_left = joined_df.column("count").map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;
    let count_right = joined_df.column("count_right").map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;

    let mut diff_series = count_left.clone() - count_right.clone();
    let diff_series = diff_series.rename("h1_h2").clone();

    let differences_df = joined_df
        .with_column(diff_series)
        .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;
    

    let sorted_df = differences_df.sort(["h1_h2"],Default::default());
    let sorted_df = sorted_df.map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;

    let upper = sorted_df
                .lazy()
                .filter(col("bin1").lt_eq(col("bin2")))
                .collect().map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;

    // Sum upper h1
    let factor_h1= upper
                        .column("count")
                        .map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?
                        .sum::<f32>()
                        .map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;
    // Sum upper h2
    let factor_h2 = upper
                        .column("count_right")
                        .map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?
                        .sum::<f32>()
                        .map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?;
    
    let factor_h1_u32 = factor_h1.ceil() as usize;
    let factor_h2_u32 = factor_h2.ceil() as usize;



    let mut diff_series = Vec::new();
    let mut bin1_stats = None;
    let mut bin2_stats = None;

    for t in 0..n_permutations{
        println!("[INFO] Starting test {}/{}", t + 1, n_permutations);
        // Perform downsample with the number of reads found in h1 and h2
        let hap1_down = downsampling(&contacts_reference,factor_h1_u32);
        let hap2_down = downsampling(&contacts_reference,factor_h2_u32);
        let hap1_down_df = array2dataframe(&hap1_down)?;
        let hap2_down_df =  array2dataframe(&hap2_down)?;
        
        //Join the two haplotype dataframes
        let stats_df = hap1_down_df.join(
            &hap2_down_df,
            ["bin1", "bin2"],
            ["bin1", "bin2"],
            JoinArgs::new(JoinType::Inner),
        )
        .map_err(|e| hdf5::Error::Internal(format!("Polars join error: {}", e)))?;
        
        // Fill the bin columns if they are empty
        if bin1_stats.is_none() || bin2_stats.is_none(){
            bin1_stats= Some(stats_df.column("bin1")
                            .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?
                            .clone());
            bin2_stats= Some(stats_df.column("bin2")
                            .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?
                            .clone());  
        }

        // Get the count columns
        let h1down_count = stats_df.column("count")
                                    .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;
        let h2down_count = stats_df.column("count_right")
                                    .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;
        
        
        let h1dc= h1down_count
                            .cast(&DataType::Float32).map_err(|e| hdf5::Error::Internal(format!("Polars cast error: {}", e)))?;
        let h2dc = h2down_count
                            .cast(&DataType::Float32).map_err(|e| hdf5::Error::Internal(format!("Polars cast error: {}", e)))?;
        
        // Calculate the difference
        let mut t_hap_diff = &h1dc - &h2dc;
        t_hap_diff.rename(&format!("trial_{}",t));
        // Add the iteration difference to the vector diff_series
        diff_series.push(t_hap_diff);
    }

    let mut all_columns = vec![
        bin1_stats.ok_or_else(|| hdf5::Error::Internal("Missing bin1".to_string()))?,
        bin2_stats.ok_or_else(|| hdf5::Error::Internal("Missing bin2".to_string()))?,
    ];
    
    all_columns.extend(diff_series
                        .into_iter()
                        .map(|s| s.clone()));
    

    let final_df = DataFrame::new(all_columns)
                            .map_err(|e| hdf5::Error::Internal(format!("Polars DataFrame error: {}", e)))?;
    

    // Merge first data frame with the final data_frame test

    let selected_columns = differences_df
            .select(["bin1", "bin2", "h1_h2"])  
            .map_err(|e| hdf5::Error::Internal(format!("Polars column addition error: {}", e)))?;

    let merged_final = selected_columns.join(
        &final_df,
        ["bin1", "bin2"],
        ["bin1", "bin2"],
        JoinArgs::new(JoinType::Inner),
    )
    .map_err(|e| hdf5::Error::Internal(format!("Polars join error: {}", e)))?;

    let original_differences: Series = merged_final.column("h1_h2")
                                            .map_err(|e| hdf5::Error::Internal(format!("Polars error: {}", e)))?
                                            .cast(&DataType::Float32)
                                            .map_err(|e| hdf5::Error::Internal(format!("Polars cast error: {}", e)))?
                                            .f32()
                                            .map_err(|e| hdf5::Error::Internal(format!("Polars conversion error: {}", e)))?
                                            .apply(|opt| opt.map(|v| v.abs())) 
                                            .into_series();

    //Vector of series where each serie represents a column 
    let test_values: Result<Vec<Series>, hdf5::Error> = merged_final
                                            .get_columns()
                                            .iter()
                                            .filter(|c| c.name().starts_with("trial_"))
                                            .map(|c| {
                                                let casted = c.cast(&DataType::Float32)
                                                    .map_err(|e| hdf5::Error::Internal(format!("Error in cast: {}", e)))?;
                                                
                                                let abs_values = casted
                                                    .f32()
                                                    .map_err(|e| hdf5::Error::Internal(format!("Error in f32: {}", e)))?
                                                    .apply(|opt| opt.map(|v| v.abs()));
                                        
                                                Ok(abs_values.into_series())
                                            })
                                            .collect();

    let test_values = test_values?;  // unwrap Result<Vec<Series>> because it was a Result<Vec<Series>, hdf5::Error

    let comparisons: Result<Vec<Series>, hdf5::Error> = test_values
                                                .iter()
                                                .map(|trial_series| {
                                                    trial_series
                                                        .gt(&original_differences)
                                                        .map(|bool_chunked| bool_chunked.into_series())
                                                        .map_err(|e| hdf5::Error::Internal(format!("Comparison error: {}", e)))
                                                })
                                                .collect();
    
    let comparisons = comparisons?;
    
    let bool_comparisons_df = DataFrame::new(comparisons)
                                                .map_err(|e| hdf5::Error::Internal(format!("polars DF creation error: {}",e)))?;
    let n_rows = bool_comparisons_df.height();

    let true_counts_vec: Vec<u32> = (0..n_rows)
                                                .map(|i| {
                                                    bool_comparisons_df
                                                        .get_columns()
                                                        .iter()
                                                        .map(|col| col.bool().unwrap().get(i).unwrap_or(false) as u32)
                                                        .sum()
                                                })
                                                .collect();
    let true_counts_series = Series::new("true_count", true_counts_vec);   
    let denom = n_permutations as f32;
    let pvalues_series = true_counts_series
            .cast(&DataType::Float32) // cast to float
            .map_err(|e| hdf5::Error::Internal(format!("Polars cast error: {}", e)))?
            .f32()
            .map_err(|e| hdf5::Error::Internal(format!("Polars conversion error: {}", e)))?
            .apply(|opt| opt.map(|v| v / denom))
            .into_series(); 

    let mut binding = pvalues_series.clone();
    let pvalues = binding.rename("p_value");
    let mut to_save = differences_df
    .hstack(&[pvalues.clone()])
    .map_err(|e| hdf5::Error::Internal(format!("Polars hstack error: {}", e)))?;
    println!("{:?}",to_save);   


    /// Save
    /// 
    let reference_path = Path::new(&args.reference);

    let name_sample = reference_path
    .file_name()                     // e.g. "Hybrid_week0_S2.merged_mq0.mcool"
    .and_then(|f| f.to_str())       // convert to &str
    .and_then(|s| s.strip_suffix(".mcool")) // remove .mcool
    .ok_or_else(|| hdf5::Error::Internal("Failed to extract sample name".to_string()))?;

    let haplodiffc_path = reference_path
    .parent()     // up 1
    .and_then(|p| p.parent())  // up 2
    .and_then(|p| p.parent())  // up 3
    .ok_or_else(|| hdf5::Error::Internal("Failed to go up 3 directories".to_string()))?
    .join("HaplodiffC").join(name_sample);

    std::fs::create_dir_all(&haplodiffc_path)
    .map_err(|e| hdf5::Error::Internal(format!("Failed to create output directory: {}", e)))?;

    let file_name = format!(
    "HaplodiffC_chr{}_chr{}_binsize_{}.csv",
        &args.chr_x, 
        &args.chr_y,
        &args.bin_size
        );

    let full_path = (haplodiffc_path).join(file_name);

    let mut file = STDFILE::create(&full_path)
    .map_err(|e| hdf5::Error::Internal(format!("Failed to create file: {}", e)))?;

    CsvWriter::new(&mut file)
    .finish(&mut to_save)
    .map_err(|e| hdf5::Error::Internal(format!("Failed to write CSV: {}", e)))?;
    println!("[INFO] File written to: {:?}", full_path);
    Ok(())

}

