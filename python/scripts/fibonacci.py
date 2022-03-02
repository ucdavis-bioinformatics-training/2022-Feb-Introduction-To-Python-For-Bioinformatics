def fib(thresh):
    n1 = 1
    n2 = 1
    retval=[1,1]

    while (n2 < thresh):
        tot = n1 + n2
        if (tot > thresh):
            break
        retval.append(tot)
        n1 = n2
        n2 = tot

    return(retval)

print(fib(1000))