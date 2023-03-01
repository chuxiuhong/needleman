use std::cmp::{max};

///
/// Return Alignment score and alignment method
#[derive(Debug)]
pub struct AlignmentResult {
    pub score: i32,
    pub res1: String,
    pub res2: String,
}

impl AlignmentResult {
    pub fn new(score: i32, res1: String, res2: String) -> AlignmentResult {
        AlignmentResult { score, res1, res2 }
    }
}

#[derive(Debug)]
pub struct KbandAlignmentResult {
    pub score: i32,
    pub res1: String,
    pub res2: String,
    pub k_band: usize,
}


impl KbandAlignmentResult {
    pub fn new(score: i32, res1: String, res2: String, k_band: usize) -> KbandAlignmentResult {
        KbandAlignmentResult{score, res1, res2, k_band}
    }
}

///
/// needleman_wunsch algorithm implement, O(n^2) complexity
///
/// s1,s2: sequence string
///
/// match_score: score which two char match
///
/// mismatch_score: score which two char mismatch
///
/// indel_score: score which a char and a gap
///
/// ignore_case: if ignore case
///
/// any_char: the char can match any char and match_socre = 0
pub fn needleman(s1: &str, s2: &str, match_score: i32, mismatch_score: i32, indel_score: i32, ignore_case: bool, any_char: char) -> AlignmentResult {
    let xlength = s1.len() + 1;
    let ylength = s2.len() + 1;
    let string1;
    let string2;
    let anychar;
    if ignore_case {
        string1 = s1.to_lowercase();
        string2 = s2.to_lowercase();
        anychar = any_char.to_lowercase().to_string();
    } else {
        string1 = s1.to_string();
        string2 = s2.to_string();
        anychar = any_char.to_string();
    };
    let s1 = string1.as_bytes();
    let s2 = string2.as_bytes();
    let mut mat = Vec::with_capacity(xlength);
    for i in 0..xlength {
        mat.push(Vec::with_capacity(ylength));
        for _j in 0..ylength {
            mat[i].push(i32::MIN);
        }
    }
    mat[0][0] = 0i32;
    for i in 1..xlength {
        mat[i][0] = i as i32 * indel_score;
    }
    for i in 1..ylength {
        mat[0][i] = i as i32 * indel_score;
    }
    for i in 1..xlength {
        for j in 1..ylength {
            mat[i][j] = if s1[i - 1] == s2[j - 1] {
                max(max(mat[i - 1][j - 1] + match_score, mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
            } else if s1[i - 1] == anychar.as_bytes()[0] || s2[j - 1] == anychar.as_bytes()[0] {
                max(max(mat[i - 1][j - 1], mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
            } else {
                max(max(mat[i - 1][j - 1] + mismatch_score, mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
            };
        }
    }
    let final_score = mat[xlength - 1][ylength - 1];
    let mut res1 = String::new();
    let mut res2 = String::new();
    let mut x = xlength - 1;
    let mut y = ylength - 1;
    while x >= 1 || y >= 1 {
        if x == 0 {
            res2.push(s2[y - 1] as char);
            res1.push('-');
            y -= 1;
        } else if y == 0 {
            res1.push(s1[x - 1] as char);
            res2.push('-');
            x -= 1;
        } else {
            let step_score = if s1[x - 1] == s2[y - 1] {
                match_score
            } else if s1[x - 1] as char == 'n' || s2[y - 1] as char == 'n' {
                0
            } else {
                mismatch_score
            };
            if mat[x - 1][y - 1] + step_score >= mat[x - 1][y] + indel_score && mat[x - 1][y - 1] + step_score > mat[x][y - 1] + indel_score {
                res1.push(s1[x - 1] as char);
                res2.push(s2[y - 1] as char);
                x -= 1;
                y -= 1;
            } else if mat[x - 1][y] > mat[x][y - 1] {
                res1.push(s1[x - 1] as char);
                res2.push('-');
                x -= 1;
            } else {
                res1.push('-');
                res2.push(s2[y - 1] as char);
                y -= 1;
            }
        }
    }
    res1 = res1.chars().rev().collect();
    res2 = res2.chars().rev().collect();
    AlignmentResult::new(final_score, res1, res2)
}


pub fn kband_needleman(s1: &str, s2: &str, match_score: i32, mismatch_score: i32, indel_score: i32, ignore_case: bool, any_char: char) -> KbandAlignmentResult {
    let xlength = s1.len() + 1;
    let ylength = s2.len() + 1;
    let string1;
    let string2;
    let anychar;
    let mut found = false;
    let mut final_score = i32::MIN;
    let mut k_band = 1;
    if ignore_case {
        string1 = s1.to_lowercase();
        string2 = s2.to_lowercase();
        anychar = any_char.to_lowercase().to_string();
    } else {
        string1 = s1.to_string();
        string2 = s2.to_string();
        anychar = any_char.to_string();
    };
    let s1 = string1.as_bytes();
    let s2 = string2.as_bytes();
    let mut mat = Vec::with_capacity(xlength);
    for i in 0..xlength {
        mat.push(Vec::with_capacity(ylength));
        for _ in 0..ylength {
            mat[i].push(i32::MIN);
        }
    }
    while !found {
        mat[0][0] = 0i32;
        for i in 0..xlength {
            for j in 0..ylength {
                mat[i][j] = i32::MIN + i32::abs(indel_score) + i32::abs(mismatch_score);
            }
        }
        for i in 1..xlength {
            mat[i][0] = i as i32 * indel_score;
        }
        for i in 1..ylength {
            mat[0][i] = i as i32 * indel_score;
        }
        for i in 1..xlength {
            //for j in max(1,i-k_band / 2)..min(ylength,i + k_band / 2 + 1) {
            for j in 1..ylength {
                if i > j {
                    if i - j > k_band / 2 { continue; }
                } else if j > i && j - i > k_band / 2 { break; }
                mat[i][j] = if s1[i - 1] == s2[j - 1] {
                    max(max(mat[i - 1][j - 1] + match_score, mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
                } else if s1[i - 1] == anychar.as_bytes()[0] || s2[j - 1] == anychar.as_bytes()[0] {
                    max(max(mat[i - 1][j - 1], mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
                } else {
                    max(max(mat[i - 1][j - 1] + mismatch_score, mat[i - 1][j] + indel_score), mat[i][j - 1] + indel_score)
                };
            }
        }
        final_score = mat[xlength - 1][ylength - 1];
        if i32::abs(final_score) <= k_band as i32 {
            found = true;
        } else {
            k_band *= 2;
        }
    }

    let mut res1 = String::new();
    let mut res2 = String::new();
    let mut x = xlength - 1;
    let mut y = ylength - 1;
    while x >= 1 || y >= 1 {
        if x == 0 {
            res2.push(s2[y - 1] as char);
            res1.push('-');
            y -= 1;
        } else if y == 0 {
            res1.push(s1[x - 1] as char);
            res2.push('-');
            x -= 1;
        } else {
            let step_score = if s1[x - 1] == s2[y - 1] {
                match_score
            } else if s1[x - 1] as char == 'n' || s2[y - 1] as char == 'n' {
                0
            } else {
                mismatch_score
            };
            if mat[x - 1][y - 1] + step_score >= mat[x - 1][y] + indel_score && mat[x - 1][y - 1] + step_score > mat[x][y - 1] + indel_score {
                res1.push(s1[x - 1] as char);
                res2.push(s2[y - 1] as char);
                x -= 1;
                y -= 1;
            } else if mat[x - 1][y] > mat[x][y - 1] {
                res1.push(s1[x - 1] as char);
                res2.push('-');
                x -= 1;
            } else {
                res1.push('-');
                res2.push(s2[y - 1] as char);
                y -= 1;
            }
        }
    }
    res1 = res1.chars().rev().collect();
    res2 = res2.chars().rev().collect();
    KbandAlignmentResult::new(final_score, res1, res2, k_band)
}