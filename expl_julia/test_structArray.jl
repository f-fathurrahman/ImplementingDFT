struct MyStruct
    a::Vector{Float64}
    b::Matrix{Float64}
end

function update_struct!(o)
    a = o.a
    b = o.b
    for i in 1:size(a,1)
        a[i] = 9*i
    end
    for j in 1:size(b,2), i in 1:size(b,1)
        b[i,j] = 10*(i + j)
    end
    return
end

function main()
    N = 5
    obj1 = MyStruct(rand(5), rand(5,3))
    
    println("Before update_struct!:")
    display(obj1.a); println()
    display(obj1.b); println()
    
    update_struct!(obj1)

    println("Before update_struct!:")
    display(obj1.a); println()
    display(obj1.b); println()
end

main()