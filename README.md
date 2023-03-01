# needleman

needleman is a needleman-wunsch algorithm implement in Rust.
## Installation

Add it to your Cargo.toml
```
needleman = "*"
```

## Usage

```rust
use needleman::needleman::needleman;
let s1 = "ACTA";
let s2 = "CGAC";
// match_score = 4, mismatch_score=-3, gap_score=-4, ignore_case=true,anychar='N'
let ag = needleman(s1,s2,4,-3,-4,true,'N');
println!("score = {}",ag.score);
println!("align1 = {}, align2 = {}",ag.res1,ag.res2);
```

Banded needleman
```rust
use needleman::needleman::kband_needleman;
let s1 = "ACTA";
let s2 = "CGAC";
let ag = kband_needleman(s1,s2,0,-1,-1,true,'N');
println!("score = {}",ag.score);
println!("align1 = {}, align2 = {}",ag.res1,ag.res2);
println!("k_band = {}",ag.k_band);
```


