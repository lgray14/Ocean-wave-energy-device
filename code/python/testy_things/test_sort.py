list1 = [3,    1,   4,   1,   5,   9,   2,   6,   5,   3,   5]
list2 = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"]

combined_lists = list(zip(list1, list2))
print(combined_lists)
sorted_lists = sorted(combined_lists, key=lambda x: x[0])
print(sorted_lists)

sorted_list2 = [item[1] for item in sorted_lists]
print(sorted_list2)