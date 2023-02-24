use std::fs;
use std::env;
use crate::needleman::needleman;
use std::time::{Duration,Instant};
mod needleman;

fn preprocess_gene(s:String) -> String{
    let mut res = String::new();
    for c in s.chars() {
        if c.is_ascii_alphabetic() {
            res.push(c);
        }
    }
    res
}




fn main() {
    let s1 = "ACTA";
    let s2 = "CGAC";
    let ag = needleman(s1,s2,4,-3,-4,true,'N');
    println!("score = {}",ag.score);
    println!("align s1 = {}, s2 = {}",ag.res1,ag.res2);
}
