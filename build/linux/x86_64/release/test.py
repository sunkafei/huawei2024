from os import system
L = [f"cases/testcase{i}.in" for i in range(0, 10, 2)] + ["lq.in"]
for filename in L:
    print(f"-------------------{filename}-------------------")
    system(f"./main {filename}")