
file = open("output_data.txt", 'r')

print(file.read())


my_list = [line.split(' , ')for line in open ("output_data.txt")]
print(my_list)

# lines = file.read().split(',')
# print(lines)
# initials = model

#
# lines = file.readlines()
#
# my_list = [line.split(' , ')for line in open ("output_data")]
# print (my_list)
#
# part1 = my_list[0:37]
# init = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000]
# final = init.extend(part1)
# print(len(final))
# print(final)

