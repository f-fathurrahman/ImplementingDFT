push!(LOAD_PATH, pwd())

using KSDFT1d

function test_isolated()
    grid = FD1dGrid( (-5.0, 5.0), 51 )
    println(grid)
end

test_isolated()