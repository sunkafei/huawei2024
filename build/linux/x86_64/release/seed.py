# baseline: 1619228
from os import listdir
max_score = {}
answer = {}
items = {}
seeds = []
for seed in listdir('data'):
    filename = f"data/{seed}"
    seed = seed.removesuffix(".csv")
    seeds.append(seed)
    with open(filename, "r") as fp:
        head = fp.readline().strip().split(',')[1:]
        score1 = fp.readline().strip().split(',')[1:]
        score2 = fp.readline().strip().split(',')[1:]
        for i in range(len(head)):
            if head[i] not in max_score:
                max_score[head[i]] = 0
                answer[head[i]] = None
            now = float(score1[i]) - float(score2[i])
            items[(head[i], seed)] = now
            if now > max_score[head[i]]:
                max_score[head[i]] = now
                answer[head[i]] = seed
score = 0
for k in max_score:
    score += max_score[k]
print(answer)
print("Score: ", score)
def calc(group):
    pos = 0
    best = 0
    for seed in seeds:
        now = 0
        for i in group:
            now += items[(i, seed)]
        if now > best:
            best = now
            pos = seed
    delta = best
    for i in group:
        delta -= items[(i, answer[i])]
    print("Best seed:", pos, sep='\t')
    print("Score:     ", best, f"   ({delta})", sep='\t')
#case1:  0<=m<50
#case2:  180<=m<200&&p[1]%2==0
#case3:  50<=m<100
#csae4:  180<=m<200&&p[1]%2==0
#case5:  180<=m<200
#case6:  150<=m<170
#case7:  350<=m<400
#csae8:  350<=m<400&&p[1]%2==0
#case9:  200<=m<350
#case10: 200<=m<350&&p[1]%2==0
calc(["case2"])