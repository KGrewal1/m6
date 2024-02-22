# M6 answers

This repo contains my answers to the M6 supervision questions.

All written questions are in the folder answers, (the raw tex and the pdf are both included)

## Code

All code is written in rust. This compiler can be installed using the instructions found here:

<https://www.rust-lang.org/tools/install>

the input files need to be pasted into a folder as `M6_files/enviar/*` : I believe this is the provided folder uncompressed and with an underscore replacing the space in the name.

## MD code

This can be run using the following

```sh
cargo r -r --bin md-code
```

this will produce both html and png plots. As png plots take around 300 ms each, this can also be run without producing them by

```sh
cargo r -r --bin md_code --no-default-features
```

the html plots are somewhat more readable, as they allow zoom: this is particularly noticeable for the O-H autocorrelation function.

### Question 4

I am not entirely sure the virial expansion has been done correctly so I have not yet done 4.3

## MC code

This can be run using the following

```sh
cargo r -r --bin md-code
```

however this is currently just a hello world shell
