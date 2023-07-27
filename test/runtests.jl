using BoxCoxTrans
using Statistics: mean, var
using Test

𝐱 = [
    185,172,424,358,157,193,439,170,155,238,194,193,356,225,378,279,835,588,280,
    206,309,444,192,205,154,168,289,246,196,645,191,613,554,429,177,205,212,166,468,
    259,393,167,534,798,172,525,2031,215,547,205,200,157,944,241,244,225,425,329,201,
    300,383,235,176,367,879,316,1082,302,204,187,307,236,292,718,354,265,213,160,597,
    157,156,226,205,362,310,404,436,204,168,358,363,209,1638,179,193,225,680,277,853,
    249,499,251,529,265,250,271,294,364,155,194,610,1271,167,173,250,185,380,183,232,
    228,230,185,263,317,204,199,578,597,493,210,162,512,509,193,395,244,569,1016,257,
    207,254,201,325,201,266,233,286,193,401,283,220,168,177,316,377,569,189,181,385,
    153,162,1099,351,263,176,300,320,271,645,237,617,188,193,212,342,330,206,169,250,
    244,150,424,357,248,1485,270,319,183,198,440,162,940,187,407,458,192,172,164,260,
    209,165,285,270,388,214,163,309,741,153,170,461,246,193,158,157,167,210,194,154,
    289,441,279,203,229,165,219,200,470,215,670,202,384,161,183,263,283,268,1383,215,
    215,279,305,432,672,286,214,1168,190,292,305,545,423,251,581,194,347,186,776,187,
    253,153,661,211,188,544,464,392,425,331,204,562,171,386,219,163,184,165,159,172,
    161,256,156,445,263,197,181,232,351,510,251,281,161,167,533,512,176,199,220,270,
    414,181,170,366,323,286,430,202,264,762,249,224,495,176,256,291,150,234,223,355,
    156,216,403,155,300,345,298,170,257,208,266,921,267,315,351,153,162,444,223,193,
    220,540,327,703,196,344,264,168,281,184,182,193,340,330,358,189,240,247,153,191,
    222,162,363,159,199,184,273,158,278,166,250,272,160,156,211,381,841,180,166,219,
    258,275,257,302,289,255,180,176,330,227,201,765,154,491,230,404,160,156,228,231,
    376,156,1195,640,607]

# precision tolerance
mytol=1e-4

# lambda tests
v = -0.991720
λ, details = BoxCoxTrans.lambda(𝐱)
@test λ ≈ v atol=mytol
λ, details = BoxCoxTrans.lambda(𝐱; method = :geomean)
@test λ ≈ v atol=mytol
λ, details = BoxCoxTrans.lambda(𝐱; method = :normal)
@test λ ≈ v atol=mytol

𝐲 = BoxCoxTrans.transform(𝐱)
@test sum(𝐲) ≈ 405.682126 atol=mytol
@test mean(𝐲) ≈ 1.0041636803675948 atol=mytol
@test var(𝐲) ≈ 2.7241853245802466e-6 atol=mytol

𝐲2 = BoxCoxTrans.transform(𝐱; scaled = true)
@test var(𝐲2) ≈ 15188.833710662804 atol=mytol

@test_throws DomainError BoxCoxTrans.transform([1,2,3,0])
@test_throws DomainError BoxCoxTrans.transform([1,2,3,-4])
@test_throws ArgumentError BoxCoxTrans.lambda(𝐱; method = :badmethod)


# BoxCoxTransformation tests
@test_throws ArgumentError BoxCoxTransformation(𝐱; method = :badmethod)

v = -0.991720
bc = BoxCoxTransformation(𝐱)
λ, _ = BoxCoxTrans.lambda(𝐱)
@test λ == bc.λ
bc = BoxCoxTransformation(𝐱; method = :geomean)
λ, _ = BoxCoxTrans.lambda(𝐱; method = :geomean)
@test λ == bc.λ
bc = BoxCoxTransformation(𝐱; method = :normal)
λ, _ = BoxCoxTrans.lambda(𝐱; method = :normal)
@test λ == bc.λ

bc = BoxCoxTransformation(𝐱)
@test BoxCoxTrans.transform(bc) == BoxCoxTrans.transform(𝐱)
@test BoxCoxTrans.transform(bc; scaled = true) == BoxCoxTrans.transform(𝐱; scaled = true)

# confidence interval
bc = BoxCoxTransformation(𝐱)
@time conf = confint(bc; alpha=0.05)
c = (-1.20795, -0.78392)
@test conf[1] ≈ c[1] atol=mytol
@test conf[2] ≈ c[2] atol=mytol

@time conf = confint(bc; alpha=0.001)
c = (-1.35975, -0.64747)
@test conf[1] ≈ c[1] atol=mytol
@test conf[2] ≈ c[2] atol=mytol
