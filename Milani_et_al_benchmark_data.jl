using MAT, Statistics, Plots, Distributions, Random, CSV, NetCDF, DataFrames,ChemStats, ProgressBars
include("Gen1D.jl")
include("GenPeak.jl")
include("GenNoise.jl")
include("GenBaseline.jl")

nominal = [0.093136565,0.000809829,3.75,1,1000,0.004]

E2 = [-0.25, 
-0.191299419,
0.004912528,
0.093136565,
0.189831804,
0.25]

S2 = [0.00048523,
0.000571715,
0.000657706,
0.000809829,
0.001035262,
0.001515892,
0.00301669]

S1 = [1.875,
3.75,
5.625,
7.5,
11.25]

PR = [0.05,
0.1,
0.5,
1,
2,
10,
20]

KT = [1, 
2,
5,
10,
100,
1000,
10000]

SW = [-0.02,
-0.01,
-0.009,
-0.008,
-0.007,
-0.006,
-0.004,
-0.002,
-0.001,
0.00,
0.001,
0.002,
0.004,
0.006,
0.007,
0.008,
0.009,
0.01,
0.02]

gridsize = 100

vars = [E2,S2,S1,PR,KT,SW]

varslb = ["E2","S2","S1","PR","KT","SW"]

Rs_lim = 3
dfpara = DataFrame(E2 = [], S2 = [], S1 = [], PR = [], KT = [],SW = [])
df4 = DataFrame(loc1D = [], loc2D =[], E2 = [], S2 = [], S1 = [], PR = [], KT = [],SW = [])
for i = 1:length(vars)
    para = vars[i]    
    temp = ones(length(para),length(vars))
    temp[:,i] = para
    filler = ones(length(para),length(vars)-1) .* nominal[1:end .!= i]'
    nm = collect(1:length(vars))[1:end .!= i]
    temp[:,nm] = filler
    append!(dfpara,DataFrame(temp,[:E2,:S2,:S1,:PR,:KT,:SW]))
end

dfaux = DataFrame(lb1D = [], hb1D = [], lb2D = [], hb2D = [])
for j = 1:size(dfpara,1)
    dt1 = 2 .* Rs_lim .* (dfpara.S1[j]+nominal[3])
    dt2 = 2 .* Rs_lim .* 60 .* (dfpara.S2[j]+nominal[2])
    l1D = collect(range(3000-dt1,3000+dt1, length=gridsize))
    l2D = collect(range(3-dt2,3+dt2, length=gridsize))
    df2 = DataFrame(loc1D = [], loc2D =[])
    sz1 = length(l1D)
    sz2 = length(l2D)
    dftm = DataFrame(loc1D = [], loc2D =[])
    for i = 1:sz2
        temp = ones(sz1,2)
        temp[:,1] = l1D
        temp[:,2] = temp[:,2] .*l2D[i]
        append!(dftm,DataFrame(temp,[:loc1D,:loc2D]))
    end
    
  
    push!(dfaux,[l1D[1],l1D[end],l2D[1],l2D[end]])
    
    temp = DataFrame(ones(length(l1D) * length(l2D),6) .* Vector(dfpara[j,:])', [:E2,:S2,:S1,:PR,:KT,:SW])
    temp2 = hcat(dftm,temp)
    df4 = vcat(df4,temp2)
end

temp = DataFrame(varval = [],varlb = [])
for i = 1:size(vars,1)
    temp2 = DataFrame(varlb = [])
    for j = 1:size(vars[i],1)
        append!(temp2,DataFrame(varlb = varslb[i]))
    end
    temp3 = hcat(DataFrame(varval = vars[i]),temp2)
    append!(temp,temp3)
end

dfaux = hcat(dfaux,temp)

## post gen mod to norm. grid size

l1D = collect(range(2910,3090,length=gridsize))
l2D = collect(range(1.62245316,4.37754684,length=gridsize))
sz1 = length(l1D)
sz2 = length(l2D)
dftm = DataFrame(loc1D = [], loc2D =[])
for j = 1:size(dfaux,1)
    for i = 1:sz2
        temp = ones(sz1,2)
        temp[:,1] = l1D
        temp[:,2] = temp[:,2] .*l2D[i]
        append!(dftm,DataFrame(temp,[:loc1D,:loc2D]))
    end
end

df4[:,1:2] = dftm 

dfaux[:,1:4] = [l1D[1] l1D[end] l2D[1] l2D[end]] .* ones(size(dfaux,1),4)


dfst = DataFrame([3000 3 nominal'], [:loc1D,:loc2D,:E2,:S2,:S1,:PR,:KT,:SW])

dfsw = DataFrame([[],[],[],[],[],[],[],[]],[:loc1D,:loc2D,:E2,:S2,:S1,:PR,:KT,:SW])
dfswm = DataFrame([[],[],[],[],[],[],[],[]],[:loc1D,:loc2D,:E2,:S2,:S1,:PR,:KT,:SW])

for i = 1:length(SW)
    append!(dfsw,DataFrame([ones(10000) .* [3000 3 nominal'][1:end-1]'  ones(10000) .* SW[i]],[:loc1D,:loc2D,:E2,:S2,:S1,:PR,:KT,:SW]))
    append!(dfswm,DataFrame([Matrix(dftm[1:10000,:]) ones(10000) .* (nominal[1:end-1])' ones(10000) .* SW[i]],[:loc1D,:loc2D,:E2,:S2,:S1,:PR,:KT,:SW]))
end

dfst2 = ones(size(df4)) .* dfst
scnbase = vcat(dfst2,dfsw)
scncomp = vcat(df4, dfswm)
auxinfo = vcat(dfaux,DataFrame([ones(size(SW)) .* Array(dfaux[1,1:4])' SW fill("SWP",length(SW))],[:lb1D,:hb1D,:lb2D,:hb2D,:varval,:varlb]))



lbls = ["Asymmetry", "Width 2D", "Width 1D", "Peak Ratio", "M factor", "Modulation shift", "parallel Modulation shift"]
nmbrs = [6,7,5,7,7,19,19]
if !isdir(pwd() * "\\simulated data")       
    mkdir(pwd() * "\\simulated data")         #check for the presence of an output dir 
    mkdir(pwd() * "\\simulated data\\$(lbls[1])")               #create output dir
    mkdir(pwd() * "\\simulated data\\$(lbls[2])")
    mkdir(pwd() * "\\simulated data\\$(lbls[3])")
    mkdir(pwd() * "\\simulated data\\$(lbls[4])")
    mkdir(pwd() * "\\simulated data\\$(lbls[5])")
    mkdir(pwd() * "\\simulated data\\$(lbls[6])")
    mkdir(pwd() * "\\simulated data\\$(lbls[7])")
    mkdir(pwd() * "\\simulated data\\Aux files")
end

CSV.write(pwd() * "\\simulated data\\Aux files\\stationarypeaks.csv",scnbase)
CSV.write(pwd() * "\\simulated data\\Aux files\\movingpeaks.csv", scncomp)
CSV.write(pwd() * "\\simulated data\\Aux files\\aux.csv", auxinfo)

Time_1D = collect(2550:7.5:3367.5)
Time_2D = collect(0:0.01:7.49)
Time_tot = collect(2550:0.01:3374.99)

orderSG = 4             #Start value for Sovitsky Golay window
ModTime = 7.5           #Modulation time
noisewindow = 100       #Window where the noise is calculated over the firs x amount of points
NppNoises = 1000           #Amount of times the std of the noisewindow is multiplied, 3 for LOD
Overlap = 80           #Minimum percentage of overlap of peaks to assign them to the same 2D peak
Fit = "ModPearsonVII" # or Fit = "Gaussian"
freq = 100
areas = zeros(Int(size(scncomp,1)),2)
chromas_2D = zeros()
pdc = zeros(20,8)
PeakDeformationMode = "ConstantArea"

catind = DataFrame(lbl =[])
for i = 1:length(nmbrs)
    append!(catind, DataFrame(lbl = fill(lbls[i],nmbrs[i].*10000)))
end 

tmp = DataFrame(val = [])
for i = 1:size(auxinfo,1)
    append!(tmp, DataFrame(val = fill(replace("$(auxinfo.varval[i])", '.' => '_'),10000)))
end  

catind = hcat(catind,tmp)
##Input parameters

for i in ProgressBar(1:length(scncomp[:,1]))

    TestTitle = "DeBug_"
    SepTech = "GC"                                                              #"GC" or "LC"
    DimensionMode = "2D"                                                        #Pick Chroma Mode by selecting "1D" or "2D"
    NumSignals = 1                                                              #Number of signals to create
    BaselineShapes_1D = 1                                                       #Simulate signal with one of the experimental baseline shapes.
    BaselineShapes_2D = 2                                                       #Simulate signal with one of the experimental baseline shapes.
    BaselineShapes_GCxGC = 6                                                    #real GCxGC baseline (only option == 6)
    BaselineDrift = 0
    Peakshapes_1D = [6, 6]                                                      #See reference table above
    Peakshapes_2D = [6, 6]                                                      #See reference table above
    Noise_sd = 0.014                                                            #Noise standaard deviatiom
    Noise_Distribution = "Normal"                                               #Noise distribution type
    Noise_Extra = 1                                                             #Noise extra
    Noise_Type = 2                                                              #Noise type, see list above
    t0_Pos_1D = 6.9                                                             #Add t0 at specified time (e.g. 3.5 min or something)
    t0_Pos_2D = 0.1                                                             #Add t0 at specified time (e.g. 3.5 min or something)
    Wdt_1D = [scnbase[i,"S1"], scncomp[i,"S1"]]                                 #half Peak width
    Wdt_2D = [scnbase[i,"S2"] .* 60, scncomp[i,"S2"] .* 60]                     #half Peak width
    m_1D = [1000,1000]                                                          #m = some number (5 for example)
    m_2D = [scnbase[i,"KT"], scncomp[i,"KT"]]                                   #m = some number (5 for example)
    As_1D =  [0,0]                                                              #As = Assymetry (extend of tailing, typical values of 0.1-0.2)
    As_2D = [scnbase[i,"E2"], scncomp[i,"E2"]]                                  #As = Assymetry (extend of tailing, typical values of 0.1-0.2)
    Cutoff = [3.5, 3.5]                                                         #Make sure this matches with the blank measurement time you have, if this exceeds the measurement time you will see the following error:
    Pos_1D = [scnbase[i,"loc1D"], scncomp[i,"loc1D"]]                           #peak position (tr)
    Pos_2D = [scnbase[i,"loc2D"], scncomp[i,"loc2D"]]                           #peak position (tr)
    PosShiftType_2D = ["Trend"]                                                 #RandomShift or Trend
    PosShift_2D = [scnbase[i,"SW"] *-1, scncomp[i,"SW"] *-1]                    #Amound of 2D shift (in )
    PosShiftLeftFrac_2D = [0.4]                                                 #fraction of left side asymmetry (0.0-1.0 , 0.5 is symmetrical)
    gD_1D = [1, 1]                                                              #gD is the Doppler (Gaussian) width
    gD_2D = [1, 1]                                                              #gD is the Doppler (Gaussian) width
    alpha_1D = [1, 1]                                                           #Alpha is the shape constant (ratio of the Lorentzian width gL to the Doppler width gD
    alpha_2D = [1, 1]                                                           #Alpha is the shape constant (ratio of the Lorentzian width gL to the Doppler width gD
    OutlierPercentage = 0                                                       #Percentage of generated peaks that will be a height significantly larger or smaller than the normal peak height based on the input data
    rng = 1337                                                                  #SEED for random generation
    Automatic = 0                                                               #-|
    Manual = 1                                                                  # |--> automatic legacy
    MinSpace = 100                                                              # |
    PeakOverlap = "Yes"    

    if PeakDeformationMode .== "ConstantHeight"
        Hgt_1D = [scnbase[i,"PR"] .* Noise_sd, scncomp[i,"PR"] .* Noise_sd] .* 2016                  #Peak Height
        Hgt_2D = [scnbase[i,"PR"] .* Noise_sd, scncomp[i,"PR"] .* Noise_sd] .* 2016 
    elseif PeakDeformationMode .== "ConstantArea"

        A = [scnbase[i,"PR"] .* Noise_sd, scncomp[i,"PR"] .* Noise_sd] .* 2016                    #Peak Height
        Hgt_1D = zeros(length(A))
        Hgt_2D = zeros(length(A))
    
        for k = 1:length(A)
            f_1D = sum(modpearson(collect(0:1:600),1,Wdt_1D[k],200,m_1D[k],As_1D[k]) .* 1) 
            f_2D = sum(modpearson(collect(0:0.01:7.5),1,Wdt_2D[k],3,m_2D[k],As_2D[k]) .* 0.01) 
            Hgt_1D[k] = (A[k] / (f_1D * f_2D))
            Hgt_2D[k] = (A[k] / (f_1D * f_2D))
        end
    else 
        error("PeakDeformationMode Unknown")
    end



    a = PosShift_2D .* ones(length(Pos_1D), length(Time_1D))
    x = Time_1D' .* ones(length(Pos_1D), length(Time_1D))
    fact = transpose(a .* (x .- Pos_1D))

    Max_Hgt_1D,
    Min_Hgt_1D,
    Max_Hgt_Outlier_1D,
    Min_Hgt_Outlier_1D,
    Max_Wdt_1D,
    Min_Wdt_1D,
    Max_m_1D,
    Min_m_1D,
    Max_As_1D,
    Min_As_1D,
    Peakdata_1D =
    calcinput(Hgt_1D, Wdt_1D, m_1D, As_1D, Pos_1D, gD_1D, alpha_1D)
    # pdc[i,:]
    chroma_2D = zeros(length(Time_1D), length(Time_2D))
    for m = 1:2
        
        THR_1D = collect(range(Time_1D[1],Time_1D[end], length= 22000))
        peaks_1D, PeaksS_1D =
        ManualGeneratePeaks(Peakshapes_1D[m], THR_1D, Peakdata_1D[:,m])
        
        temp = peaks_1D .* mean(diff(THR_1D))
        temp2 = zeros(length(Time_1D),size(peaks_1D,2))
        for q = 2:length(Time_1D)
           temp2[q,:] = sum(temp[Time_1D[q-1] .< THR_1D .< Time_1D[q],:],dims=1) 
        end
        
        peaks_1D = temp2 .* 7.5
        PeaksS_1D = sum(peaks_1D,dims = 2)


        areamod = zeros(length(PeaksS_1D),1)
        chromas = zeros(length(Time_1D), length(Time_2D))
        for k = 1:length(PeaksS_1D)
            Max_Hgt_2D,
            Min_Hgt_2D,
            Max_Hgt_Outlier_2D,
            Min_Hgt_Outlier_2D,
            Max_Wdt_2D,
            Min_Wdt_2D,
            Max_m_2D,
            Min_m_2D,
            Max_As_2D,
            Min_As_2D,
            Peakdata_2D = calcinput(
            peaks_1D[k, :],
            Wdt_2D[m],
            m_2D[m],
            As_2D[m],
            (Pos_2D[m] .+ fact[k, m]),
            gD_2D[m],
            alpha_2D[m],
            )
            peaks_2D, PeaksS_2D =
            ManualGeneratePeaks(Peakshapes_2D[m], Time_2D, Peakdata_2D)
            chromas[k, :] = PeaksS_2D
            areamod[k,1] = sum(peaks_2D./freq)
        end
            
            areas[i,m] = sum(areamod)
         
            chroma_2D = chroma_2D .+ chromas
    end
    dfs = DataFrame(Time = Time_tot, Signal = reshape(chroma_2D, :,1)[:,1])
    filename = "SimChrom_" * string(i, pad = 6) * "_" * catind.lbl[i] * "_" * catind.val[i] * ".csv"
    filepwd = pwd() * "\\simulated data\\$(catind.lbl[i])\\"
    CSV.write(filepwd * filename, dfs)
end
CSV.write(pwd() * "\\simulated data\\Aux files\\" * "areas.csv",DataFrame(areas,[:area1,:area2]))      