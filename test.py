n = 5
li = [None for i in range(n)]
def create_function(n):
    def f(x):
        return n
    return f

for i in range(n):
    li[i] = create_function(i)


print(li[2](0))