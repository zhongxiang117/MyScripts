
def find_all_rings_dfs(nbor):
    visited = [False for i in range(len(nbor))]
    dfs = []
    rings = []
    while False in visited:
        t = visited.index(False)
        visited[t] = True
        for v in nbor[t]:
            dfs.append([t,v])
        while len(dfs):
            g = dfs.pop()
            t = g[-1]
            if visited[t]:
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
            else:
                visited[t] = True
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
                    else:
                        for i in range(len(dfs)):
                            if dfs[i][-1] == t:
                                dfs[i].append(v)        # update
                        for i in range(len(g)):
                            dfs.append([*g[i:],v])      # add all
    return rings


def find_rings_dfs(nbor):
    visited = [False for i in range(len(nbor))]
    dfs = []
    rings = []
    while False in visited:
        t = visited.index(False)
        visited[t] = True
        for v in nbor[t]:
            dfs.append([t,v])
        while len(dfs):
            g = dfs.pop()
            t = g[-1]
            if not visited[t]:
                visited[t] = True
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
                    else:
                        dfs.append([*g,v])
    return rings



# nbor = {
#     0: set(),
#     1: (2,4,5),
#     2: (1,3,4),
#     3: (2,4),
#     4: (1,2,3,5),
#     5: (1,4),
# }
#
## Resembles graph like:
##     1------2
##     / \  / \
##    /___\/___\
##   5     4    3
#
## when executed:
#
# results = find_all_rings_dfs(nbor)
#
# [3, 2, 4]
# [4, 3, 2, 1]
# [5, 4, 3, 2, 1]
# [4, 2, 1]
# [5, 4, 2, 1]
# [5, 4, 1]





