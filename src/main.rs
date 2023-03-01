



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

}
