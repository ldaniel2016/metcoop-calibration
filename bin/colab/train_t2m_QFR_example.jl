#  example of training quantile function regression on MEPS data (JBB)


using Flux
using CSV
using DataFrames, Statistics
using Dates


## 
##    D A T A     O R G A N I Z A T I O N 
##

#  read data
xd = CSV.read("data/trdata_spatial_lcc_T2m_00+036.csv");

#  create new covariates
xd.dELEV = xd[:, Symbol("SG.0")] ./ 9.81  .-  xd.elev;

#  define covariates (NOTE: for simplicity . should be avoided in column names)
kens = Symbol.(:"T2.", 0:9)
ksp  = [:x, :y, :dELEV]

#  remove missing data
xd = dropmissing(xd[:, [kens; ksp; :TA; :TIME]]);

#  define logical indices for training and validation sets
tm  = Dates.unix2datetime.(Int64.(xd.TIME));
kpr = map(u -> u in [5,10,15,20,25,30], Dates.dayofmonth.(tm));
ktr = .!kpr;

#  create standardized covariates for training and validation sets (ordered members)
x = hcat( mapslices(sort, Matrix{Float32}(xd[:, kens]), dims = 2),
          Matrix{Float32}(xd[:, ksp]) );
#mean(x, dims = 1)
#std(x, dims = 1)
x = Flux.normalise(x, dims = 1); # standardization should be based on training data only!

#  observations
y = Float32.(xd.TA);



##
##    T R A I N I N G
##

#  degree of Bernstein polynomials
const degree = 8

#  quantile levels for training
const prob   = Float32.(1:10) ./ 11
const nprob  = size(prob, 1)

#  Bernstein design matrix
dbin(x, n, p) = binomial(n, x) * p^x * (1-p)^(n-x)
const B = [dbin(d, degree, p) for p in prob, d in 0:degree]   ## size: prob, degree+1

#  quantile loss functions
function qtlossBatch(b, y, B, prob)
    qt = B * b
    mean( ((y .< qt) .- prob) .* (qt .- y) ) 
end
function qtloss2(b, y, B, prob)
    qt = B * b
    mean( ((y .< qt) .- prob * ones(Float32,1,size(b,2))) .* (qt .- y) )
end
function qtlossRaw(qt, y, B, prob)
    mean( ((y .< qt) .- prob * ones(Float32,1,size(qt,2))) .* (qt .- y) )
end

#  set batch size
nb = 100

#  make training set
n      = sum(ktr)
yy     = ones(Float32, nprob, 1) * y';
probM  = prob * ones(Float32, 1, nb);
xtmp   = x[ktr, :];
ytmp   = yy[:, ktr];
trdata = [(xtmp[i, :]', ytmp[:, i]) for i in Flux.chunk(rand(1:n, n), trunc(n/nb)+1)];
trdata = trdata[1:end-1];

#  define neural network function
model  = Chain(Dense(size(x,2), 32, relu),
               Dense(32, 16, relu),
               Dense(16, degree + 1))

#  training
loss(x, y) = qtlossBatch(model(x), y, B, probM)
function evalcb()
    println("Quantile score: ", "training = ",
            round(qtloss2(model(x[ktr,:]'), yy[:,ktr], B, prob), digits = 5),
            "  validation = ",
            round(qtloss2(model(x[kpr,:]'), yy[:,kpr], B, prob), digits = 5))
end

@time Flux.@epochs 20 Flux.train!(loss, Flux.params(model), trdata, Flux.ADAM(),
                                  cb = Flux.throttle(evalcb, 10))  


#  quantile score of raw ensemble
qtRaw = mapslices(sort, Matrix(xd[:,kens])', dims = 1);
qtlossRaw(qtRaw[:,kpr], yy[:,kpr], B, prob)

#  prediction example
B * model(x[kpr,:][1:15,:]')
qtRaw[:,kpr][:,1:15]
