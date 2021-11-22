function modify01!(a)
    # create new array, a is not modified in the calling function
    for i in 1:size(a,1)
        a = a .+ 1
    end
    return
end

function modify02!(a)
    for i in 1:size(a,1)
        a .+= 1
    end
    return
end


function modify03!(a)
    for i in 1:size(a,1)
        a[i] = a[i] + size(a,1)
    end
    return
end

function main()
    a = rand(5)
    display(a); println()
    #modify01!(a)
    #modify02!(a)
    modify03!(a)
    display(a); println()
end

main()