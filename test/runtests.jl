ENV["QUIET"]= true


dir=pwd()

F0 = filter(x -> endswith(x, ".jl") && !startswith(x,"runtests"), readdir(dir))

error_files = String[]
test_times = Float64[]

for f in F0
    try
        t = @elapsed include(f)
        @info "\033[96m$f\033[0m   $t(s)"
        push!(test_times,t)
    catch
        @warn "problem with $f"
        push!(error_files,f)
    end
end

if length(error_files)>0
    @warn "problems with $error_files"
end

tt = sum(test_times)
@info "total time: $tt(s)"
