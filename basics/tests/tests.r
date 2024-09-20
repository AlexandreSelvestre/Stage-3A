setClass("Test", representation(x = "numeric", y = "numeric"))

test <- new("Test", x = 1, y = 2)
print(slotNames(test)[1])
