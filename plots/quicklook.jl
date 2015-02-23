using Gadfly

function read_datafile(filename)
    data = readdlm(filename)
    times = data[:,1]
    num_annuli = floor(length(data[1,:])/2)
    sma = data[:,2*[1:num_annuli]+1]
    masses = data[:,2*[1:num_annuli]]
    (times, sma, masses)
end

function plot_all(dirname)
    files = readdir(dirname)
    files = filter(r"^((?!flux).)*$", files)
    smas = [exp(logsma) for logsma in linspace(log(0.05), log(5), 200)]
    mmsn = [10*sma^-1.5 for sma in smas]
    layers = Dict()
    layers["MMSN"] = layer(x=smas, y=mmsn, Geom.line, color=["MMSN"])
    for file in files
        times, sma, masses = read_datafile(dirname*file)
        layers[file] = layer(x=sma[end,:], y=masses[end,:], Geom.point, color=[file])
    end
    plot(values(layers)..., Scale.x_log10, Scale.y_log10, Guide.ylabel("Surface density [g/cm^2]"), Guide.xlabel("Orbital Radius [AU]"))
end
