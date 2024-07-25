from os import system, popen
L = [f"cases/testcase{i}.in" for i in range(0, 10, 2)] + ["lq.in"]
total = 0
for filename in L:
    print(f"-------------------{filename}-------------------")
    system(f"./main {filename} 2> tmp.txt")
    with open("tmp.txt", "r") as fp:
        result = fp.readlines()
        print(*result)
        for line in result:
            if line.startswith("Score different:"):
                total += float(line.strip().split()[-1])
print("\033[1;31m" + f"Total score: {int(total)}" + "\033[0m")
# 448623