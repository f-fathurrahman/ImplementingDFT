struct MyStruct1
    Npoints::Int64
    Nspin::Int64
    data::Array{Float64}
end

struct MyStruct2
    Npoints::Int64
    Nspin::Int64    
    data::Array{Float64,2}
end

function main()
    Npoints = 5
    Nspin = 2

    in_data = rand(Npoints,Nspin)
    obj1 = MyStruct1(Npoints,Nspin,in_data)
    obj2 = MyStruct2(Npoints,Nspin,in_data)
    
    println(typeof(obj1.data))
    println(typeof(obj2.data))

    in_data2 = rand(Npoints)
    obj3 = MyStruct1(Npoints,Nspin,in_data2)
    obj4 = MyStruct2(Npoints,Nspin,reshape(in_data2,(Npoints,1))) # must use reshape

    println(typeof(obj3.data))
    println(typeof(obj4.data))
    
    println("Pass here")
end

main()