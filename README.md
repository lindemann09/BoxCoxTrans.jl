# BoxCoxTrans.jl

[![Build Status](https://travis-ci.org/tk3369/BoxCoxTrans.jl.svg?branch=master)](https://travis-ci.org/tk3369/BoxCoxTrans.jl)
[![codecov](https://codecov.io/gh/tk3369/BoxCoxTrans.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tk3369/BoxCoxTrans.jl)
[![Coverage Status](https://coveralls.io/repos/github/tk3369/BoxCoxTrans.jl/badge.svg?branch=master)](https://coveralls.io/github/tk3369/BoxCoxTrans.jl?branch=master)

This package provides an implementation of Box Cox transformation.
See [Wikipedia - Power Transform](https://en.wikipedia.org/wiki/Power_transform)
for more information.

## Installation

```
] add https://github.com/tk3369/BoxCoxTrans.jl
```

## Usage

The simplest way is to just call the `transform` function with an array of numbers.

```
julia> using Distributions, UnicodePlots, BoxCoxTrans

julia> x = rand(Gamma(2,2), 10000) .+ 1;

julia> histogram(x)
               ┌──────────────────────────────────────────┐
     (0.0,2.0] │▇▇▇▇▇▇▇▇ 852                              │
     (2.0,4.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 3554   │
     (4.0,6.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 2655            │
     (6.0,8.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1524                      │
    (8.0,10.0] │▇▇▇▇▇▇▇▇ 798                              │
   (10.0,12.0] │▇▇▇ 316                                   │
   (12.0,14.0] │▇▇ 170                                    │
   (14.0,16.0] │▇ 76                                      │
   (16.0,18.0] │ 35                                       │
   (18.0,20.0] │ 14                                       │
   (20.0,22.0] │ 5                                        │
   (22.0,24.0] │ 1                                        │
               └──────────────────────────────────────────┘

julia> histogram(BoxCoxTrans.transform(x))
             ┌────────────────────────────────────────┐
   (0.0,0.2] │▇▇ 64                                   │
   (0.2,0.4] │▇▇▇▇ 166                                │
   (0.4,0.6] │▇▇▇▇▇▇▇▇▇▇ 386                          │
   (0.6,0.8] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 559                     │
   (0.8,1.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 887             │
   (1.0,1.2] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1091      │
   (1.2,1.4] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1257  │
   (1.4,1.6] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1298 │
   (1.6,1.8] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1281 │
   (1.8,2.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1103      │
   (2.0,2.2] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 877             │
   (2.2,2.4] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 541                      │
   (2.4,2.6] │▇▇▇▇▇▇▇ 272                             │
   (2.6,2.8] │▇▇▇▇ 153                                │
   (2.8,3.0] │▇ 50                                    │
   (3.0,3.2] │ 14                                     │
   (3.2,3.4] │ 1                                      │
             └────────────────────────────────────────┘
```

You can examine the power transform parameter (λ) dervied by the program:
```
julia> BoxCoxTrans.lambda(x).value
0.013544484565969775
```

You can transfrom the data using your own λ:
```
julia> histogram(BoxCoxTrans.transform(x, 0.01))
             ┌────────────────────────────────────────┐
   (0.0,0.2] │▇▇ 64                                   │
   (0.2,0.4] │▇▇▇▇ 166                                │
   (0.4,0.6] │▇▇▇▇▇▇▇▇▇▇ 386                          │
   (0.6,0.8] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 563                     │
   (0.8,1.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 894             │
   (1.0,1.2] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1094      │
   (1.2,1.4] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1261  │
   (1.4,1.6] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1302 │
   (1.6,1.8] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1289 │
   (1.8,2.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1111      │
   (2.0,2.2] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 868             │
   (2.2,2.4] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 534                      │
   (2.4,2.6] │▇▇▇▇▇▇▇ 268                             │
   (2.6,2.8] │▇▇▇▇ 142                                │
   (2.8,3.0] │▇ 47                                    │
   (3.0,3.2] │ 11                                     │
             └────────────────────────────────────────┘
```

There's an option to scale the results by the geometric mean.
```
julia> histogram(BoxCoxTrans.transform(x; scaled = true))
               ┌────────────────────────────────────────┐
     (0.0,1.0] │▇▇ 81                                   │
     (1.0,2.0] │▇▇▇▇▇▇ 258                              │
     (2.0,3.0] │▇▇▇▇▇▇▇▇▇▇▇▇ 540                        │
     (3.0,4.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 868                 │
     (4.0,5.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1234       │
     (5.0,6.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1467  │
     (6.0,7.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1521 │
     (7.0,8.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1488  │
     (8.0,9.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 1133          │
    (9.0,10.0] │▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 806                  │
   (10.0,11.0] │▇▇▇▇▇▇▇▇ 359                            │
   (11.0,12.0] │▇▇▇▇ 188                                │
   (12.0,13.0] │▇ 50                                    │
   (13.0,14.0] │ 7                                      │
               └────────────────────────────────────────┘
```
