import random, string, sys

if len(sys.argv) < 5:
	print("provide the output file, the lenghts of the two sigmas and the lengths of the two strings, pls.")
	exit(1)

#alphabets = string.printable[:94] # till the ~ symbol
alphabets = string.ascii_lowercase

random.seed()

sigma_1 = random.sample(alphabets, int(sys.argv[2]))
sigma_2 = random.sample(alphabets, int(sys.argv[3]))

print("sample sigma 1:", sigma_1)
print("sample sigma 2:", sigma_2)

with open(sys.argv[1],'w') as f:
    f.write("".join([random.choice(sigma_1) for _ in range(int(sys.argv[4]))]))
    f.write("\n")
    f.write("".join([random.choice(sigma_2) for _ in range(int(sys.argv[5]))]))