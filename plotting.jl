using PyCall, PyPlot
using Images
using LaTeXStrings
using ThreePhotons
using Interact

"""Plots random 2 photon correlation slices"""
function plot_random_2photon_slices(c2list::Dict)
    f = figure()
    N,K,_ = Base.size(c2list[collect(keys(c2list))[1]])
    @manipulate for k1=1:K,k2=1:K,normalization=true;
        plot_single_2photon_slice(c2list, k1, k2, f, normalization)
    end
end

function plot_single_2photon_slice(c2list::Dict, k1::Int64, k2::Int64, f=figure(), normalization=true)
    withfig(f) do
        N,K,_ = Base.size(c2list[collect(keys(c2list))[1]])
        x = linspace(0,pi,N)
        for name in keys(c2list)
          slice = real(c2list[name][:,k2,k1])
          if normalization slice /= sumabs(slice) end
          plot(x,slice,"-o", label=name)
        end

        title(latexstring("2 Photon Correlation \$k_1=$k1, k_2=$k2\$"))
        xlabel(L"$\alpha$ [rad]")
        ylabel(latexstring("c_2($k1,$k2)") )
        legend()
    end
end

"""Plots random 3 photon correlation slices"""
function plot_random_3photon_slices(c3list::Dict, GaussianKernel::Float64=0.0)
    c3names = collect(keys(c3list))
    c3num = length(c3names)
    N,_,K,_,_ = Base.size(c3list[c3names[1]])
    f2 = figure()
#     x = linspace(0,pi,N)

    @manipulate for k1=1:K,k2=1:K,k3=1:K,normalization=true;
        withfig(f2) do
            axslice = subplot2grid((3,c3num),(0,0),colspan=c3num)
            aximshow = [ subplot2grid((3,c3num),(1,i-1),rowspan=2,colspan=1) for i = 1:c3num]

            for (j,name) in enumerate(keys(c3list))
                slab = real(c3list[name][:,:,k3,k2,k1])
                if normalization slab /= sumabs(slab) end

                if GaussianKernel > 0.0
                    img = convert(Images.Image,convert(Array{Float64},slab))
                    slab = imfilter(img, Kernel.gaussian(GaussianKernel));
                end

                #Plot slabs
                aximshow[j][:imshow](slab, extent=[0, pi, 0, pi], cmap="hot")
                aximshow[j][:set_xticks]([])
                aximshow[j][:get_yaxis]()[:set_visible](false)
                aximshow[j][:set_xlabel](name)
                #               axslice[:set_clim](minimum(slice), maximum(slice))
                axslice[:set_title](latexstring("\$k_1=$k1\,k_2=$k2\,k_3=$k3\$"))
                #             axslice[:set_xlabel](L"\alpha=\beta [rad]")

                axslice[:grid](false)
                axslice[:plot](diag(slab), "-o", label=name)
            end
        end
    end
end


"""Plots an array of 2D-points as 2D scatter plot"""
function plotPoints2D(points)
  scatter( [points[i][1] for i=1:length(points)], [points[i][2] for i=1:length(points)])
  # xlim(-qmax, qmax)
  # ylim(-qmax, qmax)
end

"""Plots an array of 3D-points as 3D scatter plot"""
function plotPointCloud(points)
  # points = []
  # for i = 1:4000
  #     p,rot = pointsPerOrientation(200)
  #     p = map((x)->rot*x, p)
  #     append!(points, p)
  # end

  # colors = [qmax/sumabs(points[i]) for i = 1:length(points)]
  clf()
  scatter3D( [points[i][1] for i=1:length(points)], [ points[i][2]for i=1:length(points)], [ points[i][3]for i=1:length(points)], s=0.1)
end

"""From an intensity cube, generates scattering data and plots it"""
#noise::Noise=GaussianNoise(0.0, 0.0, false)
function plot_scattering_image(intensity::CubeVolume; number_incident_photons::Integer=10000, qmax::Float64=1.0, point_size::Float64=50.0, colorfill="red", coloredge="black", background::Bool=true)

    # noisy_intensity = get_noisy_intensity(intensity, noise)
    noisy_intensity = intensity
    p,rot = pointsPerOrientation(noisy_intensity, qmax,qmax/3.0, number_incident_photons)

    #Plot underlying intensity
    if background
        r = linspace(-qmax, qmax, noisy_intensity.cubesize)
        myslice = Float64[getVolumeInterpolated(noisy_intensity, rot*[-x,y,0]) for x=r,y=r]
        myslice = (myslice).^(0.2)
        # myslice = log(max(myslice, 1e-6))
        # myslice = log(myslice)
        fig = imshow(myslice, interpolation="hermite", extent=[-qmax, qmax, -qmax, qmax], cmap="Blues")
    end

    #plot center
    fig = scatter([0.0], [0.0], marker="+", alpha=1.0, s=30.0, color="grey")
    xlim(-qmax,qmax)
    ylim(-qmax,qmax)
    gca()[:set_aspect]("equal", adjustable="box")
    #Plot scattering photons
    scatter( [p[i][2] for i=1:length(p)], [ p[i][1] for i=1:length(p)], c=colorfill, s=point_size, alpha=1.0, edgecolors=coloredge)#

    fig[:axes][:get_xaxis]()[:set_visible](false)
    fig[:axes][:get_yaxis]()[:set_visible](false)
end

function compare_c2_grid(a::C2, b::C2)
    _,_,K = Base.size(a)
    ac = complete_two_photon_correlation(a)
    bc = complete_two_photon_correlation(b)

    ratios = Float64[100.0*mean(ac[:,k2,k1]-bc[:,k2,k1]).^2/sumabs(ac[:,k2,k1]) for k1 = 1:K, k2=1:K]
    imshow(ratios, interpolation="none", extent=[1,K,K,1])
    title("Average Deviation [%]")
    colorbar()
    println("Difference: $(100.0*sumabs((ac-bc).^2)/sumabs(ac))")
end

"""Compares a list of histograms visually"""
function compare_histogram_with_theory(histograms::Dict=Dict(), N::Int64=32, K::Int64=16, normalization::Bool=false, ctype::String="c3", volumes::Dict=Dict(), L::Int64=10)
    histolist_c3 = Dict()
    histolist_c2 = Dict()
    for (name,file) in histograms
        _,c2,_,c3 = loadHistograms(K,K,file)
        histolist_c3[name] = c3
        histolist_c2[name] = c2
    end
    basis = calculate_basis(L,25, N, K)

    if ctype=="c3"
        c3list = Dict( name=>FullCorrelation_parallized(volume, basis, K, true, true) for (name,volume) in volumes)
        plot_random_3photon_slices(merge(histolist_c3,c3list); list=[random_triplet(K) for i = 1:20], normalization=normalization)
    else
        c2list = Dict( name=>twoPhotons(volume, basis, K, true, true) for (name,volume) in volumes)
        plot_random_2photon_slices(merge(histolist_c2,c2list), normalization=normalization, list=[random_doublet(K) for i = 1:20])
    end
end
