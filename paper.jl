using ThreePhotons
#Load all the code
using PyCall, PyPlot
using LaTeXStrings
@pyimport seaborn as sns

include("plotting.jl")

#Count total lines of code: git ls-files **/*.jl | xargs cat | wc -l
#Count lines for each file: git ls-files **/*.jl | xargs wc -l

#Defining estatics
sns.set_context("paper", font_scale=1.4, rc=Dict("lines.linewidth"=>2.0))
sns.set_style("white")
sns.set_style("ticks")
colormodel = (num) -> sns.cubehelix_palette(num, start=.5, rot=-.75)
colormodelarray = (num) -> [reshape(colormodel(num)[i,:],3) for i = 1:num]
colormodelarray2 = (num) -> [reshape(sns.cubehelix_palette(num)[i,:],3) for i = 1:num]

function plot_compare_intensity_cubes_plane(cube1::CubeVolume, cube2::CubeVolume;filename="paper/intensity_compare_x.pdf")
  c1 = real(max(abs(cube1.cube/sumabs(cube1.cube)), eps(Float64)))
  c2 = real(max(abs(cube2.cube/sumabs(cube2.cube)), eps(Float64)))

  cubesize = cube1.cubesize
  half = ceil(Integer, cubesize/2)
  third = ceil(Integer, cubesize/3)
  fourth = ceil(Integer, cubesize/4)
  range = half-third:half+third
  size = length(collect(range))
  newhalf = floor(Integer,size/2)

  slice1 = log(reshape(c1[range,range,half], size, size))
  slice2 = log(reshape(c2[range,range,half], size, size))

  color1 = reshape(colormodel(6)[3,:],3)
  color2 = reshape(colormodel(6)[5,:],3)

  plot_slice = function(slice)
#         fig = imshow(slice, interpolation="hermite", cmap=sns.dark_palette("purple", as_cmap=true))
      fig = imshow(log(max(exp(slice),1e-7)), interpolation="hermite", cmap="Greys")
      axis("off")
      fig[:axes][:get_xaxis]()[:set_visible](false)
      fig[:axes][:get_yaxis]()[:set_visible](false)
  end

  figure(figsize=(8,4), dpi=300)

  #First plot
  subplot(1,2,1)
  tight_layout(pad=0.4)
  plot_slice(slice1)
  axvline(newhalf, 0, newhalf, color=color1, lw=2.0, alpha=1.0)
  axhline(newhalf, 0, newhalf, color=color1, ls="--", lw=2.0, alpha=1.0)

  #Second plot
  subplot(1,2,2)
  plot_slice(slice2)
  axvline(newhalf, 0, newhalf, color=color2, lw=2.0, alpha=1.0)
  axhline(newhalf, 0, newhalf, color=color2, ls="--", lw=2.0, alpha=1.0)

  savefig(filename, bbox_inches="tight")
end

function plot_compare_intensity_cubes_linear(cube1::CubeVolume, cube2::CubeVolume; filename="paper/intensity_compare_x_slice.pdf", cube1_name="Reconstruction", cube2_name="Reference")
    c1 = real(max(abs(cube1.cube/sumabs(cube1.cube)), eps(Float64)))
    c2 = real(max(abs(cube2.cube/sumabs(cube2.cube)), eps(Float64)))

    cubesize = cube1.cubesize
    half = ceil(Integer, cubesize/2)
    third = ceil(Integer, cubesize/3)
    fourth = ceil(Integer, cubesize/4)
    range = half-third:half+third
    size = length(collect(range))
    newhalf = floor(Integer,size/2)

    slice1 = log(reshape(c1[range,range,half], size, size))
    slice2 = log(reshape(c2[range,range,half], size, size))

    color1 = reshape(colormodel(6)[3,:],3)
    color2 = reshape(colormodel(6)[5,:],3)

    #Third plot
    figure(figsize=(8,4), dpi=300)
    ylabel("log (Intensity)")
    xlabel("Wave Number k [1/Å]")
    xs = (collect(range)-half)*dr(cube1)
    xlim(minimum(xs), maximum(xs))
    ref1 = plot(xs,reshape(slice1[newhalf+1,:], size), color=color1, label="$(cube1_name) v")
    ref2 = plot(xs,reshape(slice1[:,newhalf+1], size), color=color1, ls="--", label="$(cube1_name) h")

    rec1 = plot(xs,reshape(slice2[newhalf+1,:], size), color=color2, label="$(cube2_name) v")
    rec2 = plot(xs,reshape(slice2[:,newhalf+1], size), color=color2, ls="--", label="$(cube2_name) h")

    legend(loc=1)
    savefig(filename, bbox_inches="tight")
end

"""Plots a list of intensites """
function plot_compare_intensities_linear(intensity_list::Array{Dict{String,Any},1}, filename::String)
    cube1 = intensity_list[1]["intensity"]
    cubesize = cube1.cubesize

    half = ceil(Integer, cubesize/2)
    third = ceil(Integer, cubesize/3)
    fourth = ceil(Integer, cubesize/4)
    range = half-third:half+third
    size = length(collect(range))
    newhalf = floor(Integer,size/2)

    xs = (collect(range)-half)*dr(cube1)

    plot_noise_slice = function(noise::CubeVolume, label, color, linestyle)
      cn = real(max(abs(noise.cube/sumabs(noise.cube)), eps(Float64)))
      slicen = log(reshape(cn[range,range,half], size, size))
      plot(xs,reshape(slicen[newhalf+1,:], size), color=color, ls=linestyle, label=label)
    end

    figure(figsize=(8,4), dpi=300)
    ylabel("log (Intensity)")
    xlabel("Wave Number k [1/Å]")
    xlim(minimum(xs), maximum(xs))
    # ylim(-22, maximum(cube1)+0.5)

    for entry in intensity_list
      plot_noise_slice(entry["intensity"], entry["label"], entry["color"], haskey(entry, "linestyle") ? entry["linestyle"] : "-")
    end

    legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=true, fancybox=true, framealpha=1.0)
    savefig(filename, bbox_inches="tight")
end

"""Prepares a set of results for a list of images, may include infinite result"""
function prepare_set(images::Array{Int64}, ppi::Int64, KMAX::Int64, K::Int64, N::Int64, include_infinite::Bool)
  #100p images:
  correlation_name_list = ["P$img" for img in images]
  if include_infinite push!(correlation_name_list, "L20_inf") end
  dirlist_colors = colormodelarray(length(correlation_name_list))

  correlation_path_list = merge( Dict("P$img" => "../parallel/data_generation/SH_$(ppi)p_N$(N)_K$(KMAX)_R$(KMAX).0_P$(img)/histo.dat" for img in images),
  include_infinite ? Dict("L20_inf" => "../expdata/correlations_N$(N)_K$(KMAX)_L20_inf.dat") : Dict())

  triplet_counts = Dict()
  for (name, path) in correlation_path_list
    _,_,c3ref_full,_ = loadHistograms(8,path)
    triplet_counts[name] = Dict()
    if contains(path, "inf")
        triplet_counts[name][K] = Inf
    else
        triplet_counts[name][K] = countTriplets(c3ref_full[:,:,1:K,1:K,1:K])
    end
  end

  correlation_list = Dict(name => Dict("path"=>correlation_path_list[name], "color"=>dirlist_colors[i], "linestyle"=>"-", "label_triplets"=>counts_to_string_long(triplet_counts[name][K], "triplets"), "shortlabel_triplets"=>counts_to_string_short(triplet_counts[name][K]), "label"=>counts_to_string_long(contains(name, "inf") ? Inf : ppi*name_to_image_count(name)), "shortlabel"=>counts_to_string_short(contains(name, "inf") ? Inf : ppi*name_to_image_count(name))) for (i,name) in enumerate(correlation_name_list))
  # for (i,name) in enumerate(correlation_name_list)
  #   Dict("path"=>correlation_path_list[name], "color"=>dirlist_colors[i], "linestyle"=>"-", "label_triplets"=>counts_to_string_long(triplet_counts[name][K], "triplets"), "shortlabel_triplets"=>counts_to_string_short(triplet_counts[name][K]), "label"=>counts_to_string_long(contains(name, "inf") ? Inf : ppi*name_to_image_count(name)), "shortlabel"=>counts_to_string_short(contains(name, "inf") ? Inf : ppi*name_to_image_count(name)))
  # end

  return correlation_list, correlation_name_list,triplet_counts,correlation_path_list
end

# FSC plots

function plot_sc(sc::Array{Float64}, sc_stderr::Array{Float64}, K::Integer, dq::Float64; label="", color="blue", ytext="Correlation", linestyle="-", plot_legend=true)
    if length(sc) == 0 return end
    q = get_qrange(K, dq)

    plt = plot(q, sc[1:K], linewidth=2.0, label=label, color=color, linestyle=linestyle)
    fill_between(q,sc[1:K]+sc_stderr[1:K], sc[1:K]-sc_stderr[1:K], color=PyPlot.matplotlib[:colors][:to_rgba](color, 0.5))
    xlabel("Wave Number k [1/Å]")
    ylabel(ytext)
    xlim(dq,dq*K)
    ylim(0.0, 1.0)
    if plot_legend
        l = legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=true, fancybox=true, framealpha=1.0)
    end
end

function add_resolution_axis(K::Integer, dq::Float64)
    ax = twiny()
    xlim(dq, K*dq)
    ax[:set_xticklabels](map((x)-> @sprintf("%.1f",x), 2.0*pi./ax[:get_xticks]()))
    xlabel("Radial Resolution [Å]")
end

function add_resolution_axis(qticks::Array{Float64}, K::Integer, dq::Float64)
    ax = twiny()
    xlim(dq, K*dq)
    ax[:set_xticks](qticks)
    ax[:set_xticklabels](map((x)-> @sprintf("%.1f",x), 2.0*pi./qticks))
    xlabel("Radial Resolution [Å]")
end

function plot_isc(sc_vs_triplets::Dict, filename::String, dirs::Array, dq::Float64, kcut::Int64=16)
    fig = figure(figsize=(6,5))
    plot(get_qrange(kcut, dq),0.5*ones(kcut),"g--", linewidth=2.0, color="#000000")
    for dir in dirs
        plot_sc(sc_vs_triplets[dir]["isc_nofitting"], sc_vs_triplets[dir]["isc_nofitting_err"], kcut, dq, linestyle=sc_vs_triplets[dir]["linestyle"], label=sc_vs_triplets[dir]["label"], color=sc_vs_triplets[dir]["color"], ytext="ISC")
    end
    add_resolution_axis(kcut, dq)
    tight_layout()
    if filename != ""
        savefig("$(filename)_ISC.pdf", bbox_inches="tight")
    end
    return fig
end

function plot_resolution_inset(sc_vs_triplets::Dict, fig, ax, dirs::Array, dq::Float64, kcut::Int64=16, width::String="25%", height::Float64=1.15)
    res_mean = [sc_vs_triplets[dir]["res"] for dir in dirs]
    res_err = [sc_vs_triplets[dir]["res_err"] for dir in dirs]
    res_err = [0.0 for dir in dirs] #TODO: Calculate error correctly

    inset_axes_module = pyimport("mpl_toolkits.axes_grid1.inset_locator")
    ia = inset_axes_module[:inset_axes](ax,
                        width=width, # width = 30% of parent_bbox
                        height=height, # height : 1 inch
                        loc=1)

    errorcolors = [sc_vs_triplets[dir]["color"] for dir in dirs]
    for i = 1:length(res_mean)
        errorbar(i, res_mean[i], fmt="D", yerr=res_err[i], color=errorcolors[i])
    end
    # semilogy(collect(1:length(res_mean)), res_mean, "D")
    ylim(minimum(res_mean)-1.0,maximum(res_mean)+1.3)
    xticks(collect(1:length(res_mean)), [sc_vs_triplets[dir]["shortlabel"] for dir in dirs], rotation="90", backgroundcolor=(1.0, 1.0, 1.0, 1.0), stretch="condensed", size=8.0)
    ylabel("Resolution [Å]", size=10.0)
    xlim(0,length(res_mean)+1)
#     ia[:axes][:get_xaxis]()[:set_visible](false)
end

function plot_fsc(sc_vs_triplets::Dict, filename::String, dirs::Array, dq::Float64, kcut::Int64=16, plot_inset::Bool=true)
    figure(figsize=(6,5))
    fig, ax = subplots(figsize=(6,5))
    plot(get_qrange(kcut, dq),0.5*ones(kcut),"g--", linewidth=2.0, color="#000000")
    res_mean = Dict(dir=>sc_vs_triplets[dir]["res"] for dir in dirs)

    for dir in dirs
        plot_sc(sc_vs_triplets[dir]["fsc"],sc_vs_triplets[dir]["fsc_err"], kcut, dq,linestyle=sc_vs_triplets[dir]["linestyle"], label=@sprintf("%s (%.1f Å)", sc_vs_triplets[dir]["label"], res_mean[dir]), color=sc_vs_triplets[dir]["color"], ytext="FSC")
    end
    add_resolution_axis( kcut, dq)

    if plot_inset
        plot_resolution_inset(sc_vs_triplets, fig, ax, dirs, dq, kcut)
    end

    tight_layout()

    if filename != ""
        savefig("$(filename)_FSC.pdf", bbox_inches="tight")
    end
    return fig,ax
end

function plot_isc_nofitting(sc_vs_triplets::Dict, filename::String, dirs::Array, dq::Float64, kcut::Int64=16)
    fig = figure(figsize=(6,5))
    plot(get_qrange(kcut, dq),0.5*ones(kcut),"g--", linewidth=2.0, color="#000000")
    for dir in dirs
        plot_sc(sc_vs_triplets[dir]["isc"], sc_vs_triplets[dir]["isc_err"] , kcut, dq,linestyle=sc_vs_triplets[dir]["linestyle"], label=sc_vs_triplets[dir]["label"], color=sc_vs_triplets[dir]["color"], ytext="ISC after Phasing")
    end
    add_resolution_axis(kcut, dq)
    tight_layout()
    if filename != ""
        savefig("$(filename)_ISC_after_fitting.pdf", bbox_inches="tight")
    end
    return fig
end

function plot_isc_comparision(sc_vs_triplets, filename::String, dirs::Array, dq::Float64, kcut::Int64)
    fig = figure(figsize=(10,5))
    s1 = subplot(1,2,1)
#     title("Before Phasing", y=1.2, loc="center")
    text(0.4, 0.1, "Before\nPhasing", fontsize=16, horizontalalignment="center", backgroundcolor="white")
    plot(get_qrange(kcut, dq),0.5*ones(kcut),"g--", linewidth=0.5, color="#000000")
    for dir in dirs
        plot_sc(sc_vs_triplets[dir]["isc_nofitting"], sc_vs_triplets[dir]["isc_nofitting_err"], kcut, dq, linestyle=sc_vs_triplets[dir]["linestyle"], label=sc_vs_triplets[dir]["label"], color=sc_vs_triplets[dir]["color"], ytext="ISC", plot_legend=false)
    end
    add_resolution_axis(kcut, dq)
    tight_layout()

    s2 = subplot(1,2,2)
#     title("After phasing", y=1.2, loc="center")
    text(0.4, 0.1, "After\nPhasing", fontsize=16, horizontalalignment="center", backgroundcolor="white")
    plot(get_qrange(kcut, dq),0.5*ones(kcut),"g--", linewidth=0.5, color="#000000")
    for dir in dirs
        plot_sc(sc_vs_triplets[dir]["isc"], sc_vs_triplets[dir]["isc_err"], kcut, dq,linestyle=sc_vs_triplets[dir]["linestyle"], label=sc_vs_triplets[dir]["label"], color=sc_vs_triplets[dir]["color"], ytext="", plot_legend=false)
    end
    s2[:set_yticklabels]([])
    l = legend(loc="upper center", bbox_to_anchor=(-0.05, -0.18), ncol=4, frameon=true, fancybox=true, framealpha=1.0)
    l[:get_frame]()[:set_edgecolor]("black")
    add_resolution_axis(kcut, dq)
    tight_layout()

    if filename != ""
        savefig("$(filename)_ISC_comparison.pdf", bbox_inches="tight")
    end
    return fig
end

function plot_sc_list(sc_vs_triplets::Dict, filename::String, dirs::Array, dq::Float64, kcut::Int64=16)

    #ISC plot
    fig1 = plot_isc(sc_vs_triplets, filename, dirs, dq, kcut)

    #FSC plot
    fig2,_ = plot_fsc(sc_vs_triplets, filename, dirs, dq, kcut)

    #ISC not fitting
    fig3 = plot_isc_nofitting(sc_vs_triplets, filename, dirs, dq, kcut)

    #Now alltogether
    fig3 = plot_isc_comparision(sc_vs_triplets, filename, dirs, dq, kcut)

    return fig1,fig2,fig3
end

function counts_to_string(num)
    pot = Integer(floor(log(10,num)))
    str = @sprintf("%.1f\\cdot10^{%i}",num / 10^pot, pot)
end

function counts_to_string_long(num, name="photons")
  if num == Inf
    latexstring("\\mathrm{infinite\\,$name}")
  else
    latexstring(counts_to_string(num),"\\,\\mathrm{$name}")
  end
end

function counts_to_string_short(num)
  if num == Inf
    latexstring("\\mathrm{infinite}")
  else
    latexstring(counts_to_string(num))
  end
end

function plot_crossing(kcut_vs_L::Dict, filename::String, dq::Float64, K::Int64, triplet_counts::Dict, kcut_vs_L_range=2:2:18, ppi::Int64=10)
  figure(figsize=[6, 5])

  kcut_vs_L_pictures_colors = colormodelarray(length(keys(kcut_vs_L)))
  kcut_vs_L_labels = Dict(img=>counts_to_string_long(ppi*name_to_image_count(img)) for img in keys(kcut_vs_L))
  resolutions = Dict()
  for (i,picture) in enumerate(kcut_vs_L_pictures)
      mean_resolution = [kcut_vs_L[picture][L]["res"] for L=kcut_vs_L_range]
      error_resolution = [kcut_vs_L[picture][L]["res_err"] for L=kcut_vs_L_range]
      resolutions[picture] = mean_resolution

      errorbar(collect(kcut_vs_L_range), mean_resolution, yerr=error_resolution, label=kcut_vs_L_labels[picture], color=kcut_vs_L_pictures_colors[i])
      plot(collect(kcut_vs_L_range), mean_resolution, "o", color=kcut_vs_L_pictures_colors[i])
  end

  xlabel(L"\mathrm{Expansion\ Order\ L}")
  ylabel(L"\mathrm{Resolution}\ [Å]")
  xlim(minimum(collect(kcut_vs_L_range))-1,maximum(collect(kcut_vs_L_range))+1)
  legend(loc=3)
  tight_layout()
  savefig(filename)
  return resolutions
end

function plot_fsc_and_crossing_in_one(sc_vs_triplets::Dict, filename::String, dirs::Array, dq::Float64, K::Int64, kcut_vs_L::Dict, triplet_counts::Dict, kcut_vs_L_range, ppi::Int64)

  fig, (ax1, ax2) = subplots(1,2, figsize=(12,5))

  qvec = get_qrange(K, dq)
  ax1[:plot](qvec,0.5*ones(K),"g--", linewidth=2.0, color="#000000")
  res_mean = Dict(dir=>sc_vs_triplets[dir]["res"] for dir in dirs)

  for dir in dirs
      sc,sc_stderr,linestyle,label,color = sc_vs_triplets[dir]["fsc"],sc_vs_triplets[dir]["fsc_err"],sc_vs_triplets[dir]["linestyle"],@sprintf("%s (%.1f Å)", sc_vs_triplets[dir]["label"], res_mean[dir]), sc_vs_triplets[dir]["color"]

      ax1[:plot](qvec, sc[1:K], linewidth=2.0, label=label, color=color, linestyle=linestyle)
      ax1[:fill_between](qvec,sc[1:K]+sc_stderr[1:K], sc[1:K]-sc_stderr[1:K], color=PyPlot.matplotlib[:colors][:to_rgba](color, 0.5))
  end
  ax1[:set_xlabel]("Wave Number k [1/Å]")
  ax1[:set_ylabel]("FSC")
  ax1[:set_xlim](dq,dq*K)
  ax1[:set_ylim](0.0, 1.0)

  ax3 = twiny(ax1)
  ax3[:set_xlim](dq, K*dq)
  ax3[:set_xticklabels](map((x)-> @sprintf("%.1f",x), 2.0*pi./ax3[:get_xticks]()))
  ax3[:set_xlabel]("Radial Resolution [Å]")

  plot_resolution_inset(sc_vs_triplets, fig, ax1, dirs, dq, K, "30%", 1.15)
  l = ax1[:legend](loc="upper center", bbox_to_anchor=(1.08, -0.18), ncol=4, frameon=true, fancybox=true, framealpha=1.0)
  # l[:get_frame]()[:set_facecolor]("green")
  l[:get_frame]()[:set_edgecolor]("black")

  # kcut_vs_L_labels = Dict(img=>counts_to_string_long(ppi*name_to_image_count(img)) for img in keys(kcut_vs_L))
  resolutions = Dict()
  for picture in keys(kcut_vs_L)
      mean_resolution = [kcut_vs_L[picture][L]["res"] for L=kcut_vs_L_range]
      error_resolution = [kcut_vs_L[picture][L]["res_err"] for L=kcut_vs_L_range]
      resolutions[picture] = mean_resolution

      ax2[:errorbar](collect(kcut_vs_L_range), mean_resolution, yerr=error_resolution, color=sc_vs_triplets[picture]["color"])
      ax2[:plot](collect(kcut_vs_L_range), mean_resolution, "o", color=sc_vs_triplets[picture]["color"])
  end

  ax2[:set_xlabel](L"\mathrm{Expansion\ Order\ L}")
  ax2[:set_ylabel](L"\mathrm{Resolution}\ [Å]")
  ax2[:set_xlim](minimum(collect(kcut_vs_L_range))-1,maximum(collect(kcut_vs_L_range))+1)
  tight_layout()
  savefig(filename, bbox_inches="tight")
end

#-------------------------------------------------------------

function plot_noise_resolution(fig, noise_mean_crossing, noise_error_crossing, gammavals, noise_colors, resolution_ylim_spacing=1.0)
#     errorbar(collect(gammarange), noise_mean_crossing, yerr=noise_error_crossing, color=noise_colors)
#     scatter(collect(gammarange), noise_mean_crossing, marker="o", color=noise_colors)
  for i = 1:length(gammavals)
      errorbar(gammavals[i], noise_mean_crossing[i], fmt="D", yerr=noise_error_crossing[i], color=noise_colors[i])
  end
  ylabel(L"\mathrm{Resolution\ [Å]}")
  xlabel(L"\mathrm{Noise\ Level}\ \gamma", backgroundcolor="white")
  xlim(minimum(gammavals)-0.1,maximum(gammavals)+0.1)

  # patches = pyimport("matplotlib.patches")
  # ax[:add_patch](patches[:Rectangle]((0.25, 3.0), 0.25, 7.0, facecolor="grey"))#alpha=0.5

  fig[:set_xticks](gammavals)
  fig[:set_xticklabels](gammavals)
  ylim(minimum(noise_mean_crossing)-resolution_ylim_spacing, maximum(noise_mean_crossing)+resolution_ylim_spacing)
end

function plot_noise_sc(noise_sc, K::Int64, dq::Float64, filename::String="paper/noise_summary.pdf", gammavals::Vector{Float64}=[0.0, 0.25, 0.5])
  noise_mean_crossing = [noise_sc[gamma]["res"] for gamma in gammavals]
  noise_error_crossing = [noise_sc[gamma]["res_err"] for gamma in gammavals]
  noise_colors = reverse(colormodelarray(length(noise_sc)))

  for (i,gamma) in enumerate(gammavals)
      noise_sc[gamma]["linestyle"] = "-"
      noise_sc[gamma]["label"] = latexstring("\$\\gamma=$gamma\$")
      noise_sc[gamma]["color"] = noise_colors[i]
  end

  # fig, ax = subplots(figsize=[6, 5])
  #plot_fsc(sc_vs_triplets, filename::String, dirs, dq::Float64, kcut::Int64=16, plot_inset::Bool=true)
  fig, ax = plot_fsc(noise_sc, "", gammavals, dq, K+10, false)
  inset_axes_module = pyimport("mpl_toolkits.axes_grid1.inset_locator")
  ia = inset_axes_module[:inset_axes](ax,
                          width="35%", # width = 30% of parent_bbox
                          height=1.30, # height : 1 inch
                          loc=1)
  plot_noise_resolution(ia, noise_mean_crossing, noise_error_crossing, gammavals, noise_colors)
  savefig(filename, bbox_inches="tight")
end

"""Plots the resolution of noisy runs for various sigmas as a function of gamma (noise level)"""
function plot_noise_various_sigma(noise_sigmas, sigmavals::Vector{Float64}, K::Int64, dq::Float64, gammavals::Vector{Float64}, filename::String)
  fig = figure(figsize=[6.4, 5])
  resolution_ylim_spacing = 0.1
  res_min = 20.0
  res_max = 0.0
  for (i,sigma) in enumerate(sigmavals)
    noise_sc = noise_sigmas[sigma]
    noise_mean_crossing = [noise_sc[gamma]["res"] for gamma in gammavals]
    noise_error_crossing = [noise_sc[gamma]["res_err"] for gamma in gammavals]
    # noise_colors = reverse(sns.color_palette(sigma_color[sigma],length(noise_sc)))

    errorbar(gammavals, noise_mean_crossing, fmt="D-", yerr=noise_error_crossing, label=latexstring("\\sigma=$sigma"), color=colormodelarray2(length(sigmavals))[i])#, color=noise_colors)

    # fig[:set_xticks](gammavals)
    # fig[:set_xticklabels](gammavals)
    res_min = minimum(noise_mean_crossing)-resolution_ylim_spacing < res_min ? minimum(noise_mean_crossing)-resolution_ylim_spacing : res_min
    res_max = maximum(noise_mean_crossing)+resolution_ylim_spacing > res_max ? maximum(noise_mean_crossing)+resolution_ylim_spacing : res_max
  end
  ylabel(L"\mathrm{Resolution\ [Å]}")
  xlabel(L"\mathrm{Noise\ Level}\ \gamma\ [\%]", backgroundcolor="white")
  ylim(res_min-0.5, res_max+0.5)
  xlim(minimum(gammavals)-0.025,maximum(gammavals)+0.025)
  xticks(gammavals)
  ax = gca()
  ax[:set_xticklabels](map((x)-> @sprintf("%d",100.0*x), gammavals))
  l = legend(loc=2)
  l[:get_frame]()[:set_edgecolor]("black")
  savefig(filename, bbox_inches="tight")
end

function calculate_ISC_vs_L(intensity::SphericalHarmonicsVolume, num_tries::Int64 = 10, ISC_vs_L_K::Int64 = 35, ISC_vs_L_range::Range = 8:2:20)
    model_isc = Dict()

    for l = ISC_vs_L_range
        model_isc[l] = Dict()
        isc = Any[]
        reduced_resolution_intensity = deleteTerms(intensity, intensity.kmax, l)
        for i = 1:num_tries
            phi,theta,gamma = get_euler_angles(random_rotation(3))
            rcoeff_A = rotateStructure(reduced_resolution_intensity, theta, phi, gamma, intensity.kmax, l)
            rcoeff_B = rotateStructure(intensity, theta, phi, gamma, intensity.kmax, intensity.lmax-1)
    #         push!(isc, shell_correlation_ISC(rcoeff_A, rcoeff_B))
            push!(isc, shell_correlation_ISC(getCube(rcoeff_A), getCube(rcoeff_B)))
        end
        model_isc[l]["all"] = isc
        model_isc[l]["mean"] = [mean(Float64[isc[i][k] for i=1:num_tries]) for k=1:intensity.kmax]#mean(isc)
        model_isc[l]["std"] = [std(Float64[isc[i][k] for i=1:num_tries]) for k=1:intensity.kmax]
        println("Processed $l")
    end
    return ISC_vs_L_K,ISC_vs_L_range,model_isc
end

function plot_ISC_vs_L(ISC_vs_L_K::Int64, ISC_vs_L_range::Range, model_isc::Dict, dq::Float64)
    figure(figsize=(6,5))
    for l = ISC_vs_L_range
        errorbar(get_qrange(ISC_vs_L_K, dq),model_isc[l]["mean"][1:ISC_vs_L_K],yerr=model_isc[l]["std"][1:ISC_vs_L_K], label=latexstring("\$L=$l\$"), color=colormodelarray(8)[round(Integer, (l-(minimum(ISC_vs_L_range)-2))/2)])
    end
    xlabel("Wave Number k [1/Å]")
    ylabel("ISC")
    # title("ISC between Full-Res Intensity and Model-based Intensity")
    legend(loc=3, fontsize=11)
    # ylim(0.2 ,1.05)
    xlim(dq,ISC_vs_L_K*dq)
    tight_layout()
    savefig("paper/ISC_vs_L.pdf", bbox_inches="tight")
end

function plot_optimal_parameter(parameter_list::Array, KMAX::Int64)
    figure(figsize=(6,5))
    Krange = 20:2:min(KMAX,40)
    Lrange = 12:2:floor(Int64,(KMAX-1)/2)
    colors = colormodelarray(length(collect(Krange)))#[Integer((K-10)/2+1)]
    #Calculating
    for (i,K) in enumerate(Krange)

        plot(collect(Lrange), Float64[get_optimal(parameter_list,K, L)["res_optimal"] for L in Lrange], label="K = $K", color=colors[i])
        y = get_optimal(parameter_list,K, maximum(Lrange))["res_optimal"]
        if i < 6
            text(maximum(Lrange)+0.25, y, "K=$K", fontsize=12, backgroundcolor="white")
        end
    end
    title(latexstring("K_{\\rm max}=$KMAX"))
    xlim(minimum(Lrange),maximum(Lrange))
    xticks(collect(Lrange))
    xlabel(L"\mathrm{Expansion\ Order\ L}")
    ylabel(L"\mathrm{Resolution\ [Å]}")
    legend(ncol=2)
    savefig("paper/res_vs_L_various_K_KMAX$(KMAX).pdf", bbox_inches="tight")
end

#Cross sections (pretty sure this is the right number - cross checked by various sites)
#@ 5keV (source: http://physics.nist.gov/cgi-bin/Xcom/xcom3_1 or http://aip.scitation.org/doi/pdf/10.1063/1.555523 TABLE 2)
barn_to_nm2 = (x) -> x* 1e-24 * 1.0e14 # [barn/atom] to [cm^2/atom] and then [cm^2/atom] to [nm^2/atom]
cross_sections_5kev = Dict("H"=>barn_to_nm2(0.13), "C"=>barn_to_nm2(7.21), "N"=>barn_to_nm2(11.1), "O"=>barn_to_nm2(16.4))
cross_sections_incoh_5kev = Dict("H"=>barn_to_nm2(0.51), "C"=>barn_to_nm2(1.98), "N"=>barn_to_nm2(2.23), "O"=>barn_to_nm2(2.32))

#@ 10 keV
# cross_sections_10kev = Dict("C"=>barn_to_nm2(3.248), "N"=>barn_to_nm2(4.722), "O"=>barn_to_nm2(6.805))

#Parameters taken from: A comprehensive simulation framework for imaging single particles and biomolecules at the European X-ray Free-Electron Laser, http://www.nature.com/articles/srep24791
"""
area: photons per nm^2 area
"""
function calculate_scattered_photons(pdb::String="../structures/crambin.pdb"; area::Float64=(100.0/2.0)^2*pi, beam_photons::Float64=5.0e11)
    fluence_CFEL = beam_photons / area
    println("Fluence ($(fluence_CFEL)), area ($(area)), beam photons ($beam_photons)")
    println("Processing $pdb")
    #Atom numbers
    molecule = loadPDB(pdb)

    count = (t) -> sum([molecule[i].t==t for i = 1:length(molecule)])
    atom_numbers = Dict( "C"=>count("C"), "N"=>count("N"), "O"=>count("O"))

    scattered_photons_coherent = fluence_CFEL * sum(atom_numbers[t] * cross_sections_5kev[t] for t in keys(atom_numbers) )
    println("Coherent photons: $scattered_photons_coherent")
    scattered_photons_incoherent = fluence_CFEL * sum(atom_numbers[t] * cross_sections_incoh_5kev[t] for t in keys(atom_numbers) )
    println("Incoherent photons: $scattered_photons_incoherent")
end

function calculate_scattered_photons_water(numer_water_atoms::Int64; area::Float64=(100.0/2.0)^2*pi, beam_photons::Float64=5.0e11)
    fluence_CFEL = beam_photons / area
    println("Fluence ($(fluence_CFEL)), area ($(area)), beam photons ($beam_photons)")
    atom_numbers = Dict( "C"=>0, "N"=>0, "O"=>numer_water_atoms)

    scattered_photons_coherent = fluence_CFEL * sum(atom_numbers[t] * cross_sections_5kev[t] for t in keys(atom_numbers) )
    println("Coherent photons: $scattered_photons_coherent")
    scattered_photons_incoherent = fluence_CFEL * sum(atom_numbers[t] * cross_sections_incoh_5kev[t] for t in keys(atom_numbers) )
    println("Incoherent photons: $scattered_photons_incoherent")
end

function name_to_image_count(name)
  return parse(Int64,replace(name,"P",""))
end

"""Returns a list of images and corresponding triplet counts"""
function images_and_triplet_counts(correlation_name_list, triplet_counts, K::Int64=26)
  images = Int64[]
  triplets = Int64[]
  for name in correlation_name_list
      if contains(name, "inf") == false
          push!(images, parse(Int64,replace(name,"P","")))
          push!(triplets, round(triplet_counts[name][K]))
      end
  end
  return images,triplets
end

function plot_data_acquiration(images_10p::Array{Int64}, triplet_counts_10p::Array{Int64})

  # figure(figsize=[6.5, 6])

  f, (ax1, ax2) = subplots(2, sharex=true, figsize=[6.5, 5])
  # subplot(2,1,1)

  data_acquisition_colormodel = (num) -> sns.color_palette("Greys",num)

  ax1[:loglog](images_10p, triplet_counts_10p, "-o", label="10 Photons", color=data_acquisition_colormodel(4)[1])
  ax1[:loglog](images_10p, [calculate_expected_triplets(img, 25) for img in images_10p], "-o", label="25 Photons", color=data_acquisition_colormodel(4)[2])
  ax1[:loglog](images_10p, [calculate_expected_triplets(img, 50) for img in images_10p], "-o", label="50 Photons", color=data_acquisition_colormodel(4)[3])
  ax1[:loglog](images_10p, [calculate_expected_triplets(img, 100) for img in images_10p], "-o", label="100 Photons", color=data_acquisition_colormodel(4)[4])
  ax1[:set_ylabel]("# Triplets")
  ax1[:tick_params](which="both", width=1.25)


  ax_top = ax1[:twiny]()
  ax_top[:set_xlabel]("Data Acquisition Time [min]")
  ax_top[:set_xscale]("log")
  ax_top[:loglog](images_10p, triplet_counts_10p, "-o", color=data_acquisition_colormodel(4)[1])
  ax_top[:set_xticks](images_10p)
  ax_top[:set_xticklabels](map((x)-> @sprintf("%1.0f",minutues_to_measure(x, 27000, 0.1)), ax_top[:get_xticks]()))
  # subplot(2,1,2)

  ax2[:loglog](images_10p, 10*images_10p, "-o", label="10 Photons", color=data_acquisition_colormodel(4)[1])
  ax2[:loglog](images_10p, 25*images_10p, "-o", label="25 Photons", color=data_acquisition_colormodel(4)[2])
  ax2[:loglog](images_10p, 50*images_10p, "-o", label="50 Photons", color=data_acquisition_colormodel(4)[3])
  ax2[:loglog](images_10p, 100*images_10p, "-o", label="100 Photons", color=data_acquisition_colormodel(4)[4])
  ax2[:set_xlabel]("# Images")
  ax2[:set_ylabel]("# Photons")
  ax2[:tick_params](which="both", width=1.25)

  # ax2_top = ax2[:twiny]()
  # # ax2_top[:set_xlabel]("Data Acquisition Time [min]")
  # ax2_top[:set_xscale]("log")
  # ax2_top[:loglog](images_10p, triplet_counts_10p, "-o", color=data_acquisition_colormodel(4)[1])
  # ax2_top[:set_xticks](images_10p)
  # ax2_top[:set_xticklabels](map((x)-> @sprintf("%1.0f",minutues_to_measure(x, 27000, 0.1)), ax_top[:get_xticks]()))

  ax2[:legend](loc=4, fontsize=11)

  f[:subplots_adjust](hspace=0.1)

  savefig("paper/data_acquisition.pdf")
end

function calculate_optimal_parameter_ratio(intensity::CubeVolume, unknown_range::Range=1000:500:7000)
  ratios = Float64[]
  for num_unknown=unknown_range
    bestr = 0
    bestc = 0
    bestK = 0
    bestL = 0
    LMAX = Integer(KMAX/2)
    for L = 10:2:24
        num_coefficients = 0.5*L^2+1.5L+1
        K = round(Int64, num_unknown / num_coefficients)

        intensity_quality = cubeToSphericalHarmonics(intensity, K, Int64(L+1))
        intensity_quality_cube = getCube(intensity_quality, intensity.cubesize)
        # rmax = K*(2*pi/intensity.rmax)
        c = sumabs(shell_correlation_ISC(intensity, intensity_quality_cube))
        # println("L=$L, K=$K, dq=$(dq(intensity_quality)), c=$(c)")
        if c > bestc
            bestc = c
            bestr = K/L
            bestK = K
            bestL = L
        end
    end
    push!(ratios,bestr)
    println("bestL: $bestL, bestK: $bestK, Unknowns: $num_unknown - Bestratio: $bestr ($bestc)")
  end

  num_unknowns = [collect(unknown_range)]
  axhline(mean(ratios), color=colormodelarray(6)[3])
  scatter(num_unknowns, ratios, s=60.0, marker="D", color=colormodelarray(6)[5])
  mr = @sprintf("%.2f",mean(ratios))
  text(minimum(unknown_range)+50, mean(ratios)-0.3, "mean ratio = $(mr)", color=colormodelarray(6)[3])
  xlim(minimum(unknown_range)-100, maximum(unknown_range)+100)
  xlabel("number of unknowns")
  ylabel("ratio K/L")
  ylim(0.0, 2.0)
  tight_layout()
  savefig("paper/ratio_vs_unknowns.pdf", bbox_inches="tight")
end


function load_linear_independant_triplets_data(ppi_list::Array{Int64}, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume, K::Int64=26, L::Int64=18, KMAX::Int64=38, N::Int64=32)
  # triplet_counts::Dict,sc_vs_triplets::Dict, correlation_name_list::Dict
  triplet_counts_ppi = Dict()
  correlation_name_list_ppi = Dict()
  sc_vs_triplets_ppi = Dict()

  for ppi in ppi_list
    println("Loading $ppi photons.")
    correlation_list,correlation_name_list,triplet_counts = prepare_set(calculate_images_ppi(ppi), ppi, KMAX, K, N, false)
    triplet_counts_ppi[ppi] = triplet_counts
    correlation_name_list_ppi[ppi] = correlation_name_list

    sc_vs_triplets_ppi[ppi] = Dict( name => merge(info, load_runs("../parallel/paper_res_vs_pictures_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_L$(L)_0.99998/$(name)", densityCube, fourierCube, intensityCube; range=1000:1019)) for (name, info) in correlation_list)
  end
  return triplet_counts_ppi, correlation_name_list_ppi, sc_vs_triplets_ppi
end

function plot_linear_independant_triplets(ppi_list::Array{Int64}, triplet_counts_ppi::Dict, correlation_name_list_ppi::Dict, sc_vs_triplets_ppi::Dict, K::Int64=26)
  fig = figure(figsize=(6,4))
  # ax1 = gca()[:xaxis]()
  # ax = subplot(1,2,1)
  xscale("log", nonposx="clip")
  # ax[:set_yscale]("log", nonposy="clip")
  for (i,ppi) in enumerate(ppi_list)
    errorbar( [ppi*name_to_image_count(dic) for dic in correlation_name_list_ppi[ppi]],
              [sc_vs_triplets_ppi[ppi][dic]["res"] for dic in correlation_name_list_ppi[ppi]],
              yerr=[sc_vs_triplets_ppi[ppi][dic]["res_err"] for dic in correlation_name_list_ppi[ppi]],
              fmt="D-", label="$(ppi)p", color=colormodelarray(length(ppi_list))[i])
  end
  legend()
  xlabel("Total Number of Photons")
  ylabel("Resolution [Å]")

  # colors = [reshape(sns.cubehelix_palette(3)[i,:],3) for i = 1:3] #sns.color_palette("Paired")[1:2:5]
  # resolutions_list = []
  # triplets_list = []
  # considered_points = [Dict(10=>3,25=>3,50=>4,100=>4),Dict(10=>4,25=>5,50=>5, 100=>5), Dict(10=>6,25=>6,50=>6,100=>6)]
  # for (i,points) in enumerate(considered_points)
  #   resolutions = Float64[]
  #   triplets = Float64[]
  #   for ppi in sort(collect(keys(points)))
  #     name = correlation_name_list_ppi[ppi][points[ppi]]
  #     push!(resolutions, sc_vs_triplets_ppi[ppi][name]["res"])
  #     push!(triplets, triplet_counts_ppi[ppi][name][K])
  #   end
  #   push!(resolutions_list, resolutions)
  #   push!(triplets_list, triplets)
  #   axhline(mean(resolutions), color=colors[i])
  # end

  # subplot(1,2,2)
  # for (i,triplets) in enumerate(triplets_list)
  #   ppi_x = sort(collect(keys(considered_points[i])))
  #   semilogy(ppi_x, triplets,"D-", label=@sprintf("Resolution = %1.1f A",mean(resolutions_list[i])), color=colors[i])
  # end
  # legend()
  # xlabel("Photons per Image (PPI)")
  # ylabel("Required Triplets")
  savefig("paper/resolution_vs_photons_various_ppi.pdf", bbox_inches="tight")
end

function plot_resolution_vs_ppi_various_total_photons(ppi_list::Array{Int64}, triplet_counts_ppi::Dict, correlation_name_list_ppi::Dict, sc_vs_triplets_ppi::Dict, K::Int64=26)
  num_photons_ppi = Dict(ppi=>calculate_images_ppi(ppi)*ppi for ppi in ppi_list)
  base_ppi = ppi_list[1]
  for (i,num_photons) in enumerate(num_photons_ppi[base_ppi][1:3])
    ppis = Float64[]
    resolutions = Float64[]
    differences = Float64[]
    push!(ppis, base_ppi)
    push!(resolutions, sc_vs_triplets_ppi[base_ppi]["P$(Integer(num_photons/base_ppi))"]["res"])
    for ppi in ppi_list[2:end]
      similar_photons = sort(num_photons_ppi[ppi], lt=(a,b)-> abs(a-num_photons) < abs(b-num_photons))[1]
      similar_photons_i = findin(num_photons_ppi[ppi],similar_photons)
      push!(ppis,ppi)
      push!(resolutions,sc_vs_triplets_ppi[ppi]["P$(Integer(similar_photons/ppi))"]["res"])
      push!(differences, abs(similar_photons-num_photons))
    end
    plot(ppis, resolutions, label=latexstring("\\mathrm{approx.}\\,$(counts_to_string(num_photons))\\pm $(counts_to_string(maximum(differences)))\\,\\mathrm{photons}"), "D-", color=colormodelarray2(3)[i])
  end
  xticks(ppi_list)
  xlabel("Photons per Image")
  ylabel("Achieved Resolution [Å]")
  legend(loc=1)
  savefig("paper/resolution_vs_ppi.pdf", bbox_inches="tight")
end
