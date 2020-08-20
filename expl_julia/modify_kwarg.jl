function myfunc!(x; y=1.1)
    y = 2*y
    println("In myfunc! y = ", y)
    z = x + y + 1
    return z
end

function myfunc2!(x; y=[1.1])
    y[1] = 2*y[1]
    println("In myfunc2! y = ", y)
    z = x + y[1] + 1
    return z
end

function main()
    x = 2.0
    y = 3.1

    z = myfunc!(x, y=y)
    println("y = ", y)
    println("z = ", z)

    y = 3.1
    z = myfunc2!(x, y=[y])
    println("y = ", y)
    println("z = ", z)

    y = [3.1]
    z = myfunc2!(x, y=y)
    println("y = ", y)
    println("z = ", z)
end

main()