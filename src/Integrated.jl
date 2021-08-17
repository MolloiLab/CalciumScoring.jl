abstract type CalciumScore end

struct IntegratedNumber{T} <: CalciumScore
    vol
    I
    N
    S_Bkg
    S_Obj
end

function IntegratedNumber(
        vol,
        I=sum(vol), 
        N=length(vol), 
        S_Bkg, 
        S_Obj
        )

    IntegratedNumber{typeof(N)}(I, N, S_Bkg, S_Obj)
end

struct IntegratedVolume{T} <: CalciumScore
    vol
    I
    N
    S_Bkg
    S_Obj
    size
end

function IntegratedVolume(
    vol,
    I=sum(vol), 
    N=length(vol), 
    S_Bkg, 
    S_Obj,
    size
    )

    IntegratedVolume{typeof(N)}(I, N, S_Bkg, S_Obj, size)
end

struct IntegratedMass{T} <: CalciumScore
    vol
    I
    N
    S_Bkg
    S_Obj
    size
    ρ
end

function IntegratedMass(
    vol,
    I=sum(vol), 
    N=length(vol), 
    S_Bkg, 
    S_Obj,
    size,
    ρ
    )

    IntegratedMass{typeof(N)}(I, N, S_Bkg, S_Obj, size, ρ)
end

function score(vol, algorithm::IntegratedNumber)
    I = algorithm.I
    N = algorithm.N
    S_Bkg = algorithm.S_Bkg
    S_Obj = algorithm.S_Obj

    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    return N_Obj
end

function score(vol, algorithm::IntegratedVolume)
    I = algorithm.I
    N = algorithm.N
    S_Bkg = algorithm.S_Bkg
    S_Obj = algorithm.S_Obj
    size = algorithm.size
    
    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    return Vol_Obj
end

function score(vol, algorithm::IntegratedMass)
    I = algorithm.I
    N = algorithm.N
    S_Bkg = algorithm.S_Bkg
    S_Obj = algorithm.S_Obj
    size = algorithm.size
    ρ = algorithm.ρ
    
    N_Obj = (I - (N * S_Bkg)) / (S_Obj - S_Bkg)
    Vol_Obj = N_Obj * size[1] * size[2] * size[3]
    Mass_Obj = Vol_Obj * ρ
    return Mass_Obj
end 
