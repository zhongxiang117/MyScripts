

def is_same_ring(al,bl):
    if not isinstance(al,list) or not isinstance(bl,list) or len(al) != len(bl):
        return False
    if len(set(al)) != len(al) or len(set(bl)) != len(bl):
        print('Warning: wrong input: duplicate')
        return False
    if len(al) <= 2:
        print('Warning: wrong input: length is 2')
        return False
    n = len(al)
    for i in range(n-1):
        v = al[i]
        if v not in bl:
            return False
        j = bl.index(v)
        bw = bl[j-1]
        if j >= n-1:
            be = bl[0]
        else:
            be = bl[j+1]
        aw = al[i+1]
        ae = al[i-1]
        if not ((aw==bw and ae==be) or (aw==be and ae==bw)):
            return False
    return True


def find_all_rings_dfs(mbors):
    """
    format:
        1--2----3
           |    |
           5----4
        nbors = {0:[], 1:[2,], 2:[3,5], 3:[2,4], 4:[3,5], 5:[2,4]}
        
        `0`: as the starting engine
    """
    # clean up
    nbors = {}
    for k,v in mbors.items():
        g = list(set(v))
        if k in g: g.remove(k)
        nbors[k] = g
    if 0 not in nbors: nbors[0] = []

    visited = [False for i in range(len(nbors))]
    dfs = []
    rings = []
    done = []
    while False in visited:
        t = visited.index(False)
        visited[t] = True
        for v in nbors[t]:
            dfs.append([t,v])
        while len(dfs):
            g = dfs.pop()
            if g in done:
                continue
            done.append(g)
            t = g[-1]
            if visited[t]:
                if g[-1] in g[:-1]:
                    i = g.index(g[-1])
                    u = g[i+1:]
                    bo = True
                    for r in rings:
                        if is_same_ring(r,u):
                            bo = False
                            break
                    if bo:
                        rings.append(u)
                    g = g[i+2:]         # update
                    dfs.append(g)
                for v in nbors[g[-1]]:
                    if v in g:
                        if v == g[-2]: continue
                        i = g.index(v)
                        u = g[i:]
                        bo = True
                        for r in rings:
                            if is_same_ring(r,u):
                                bo = False
                                break
                        if bo:
                            rings.append(u)
                        if len(g)-i > 3:
                            dfs.append([*g[i+2:],v])    # care `+2`
                    else:
                        if len(nbors[v]) != 1:
                            dfs.append([*g,v])
            else:
                visited[t] = True
                for v in nbors[t]:
                    for p in dfs:           # internal update
                        if p[-1] == t and p[-2] != v:
                            p.append(v)
                    if g[-2] != v:
                        dfs.append([*g,v])  # outside update
    return rings


# 3--2-1---4-----9
#    |/   /  \
#    5---8-7--6
nbors = {
    0:[], 1:[2,4,5], 2:[1,3,5], 3:[2,], 4:[1,6,8,9], 5:[1,2,8], 6:[4,7],
    7:[6,8], 8:[4,5,7], 9:[4,]
}
print(find_all_rings_dfs(nbors))

# 1--2--3   4
#  \  \/
#   6--5--7
nbors = {0:[], 1:[2,6], 2:[1,3,5], 3:[2,5], 4:[], 5:[2,3,6,7], 6:[1,5], 7:[5,]}
print(find_all_rings_dfs(nbors))

# 3--2-1---4-----9
#    |/   /  \
#    5---8-7--6
nbors = {
    0:[], 1:[2,4,5], 2:[1,3,5], 3:[2,], 4:[1,6,8,9], 5:[1,2,8], 6:[4,7],
    7:[6,8], 8:[4,5,7], 9:[4,]
}
print(find_all_rings_dfs(nbors))

# 6--9-5---4-----7
#    |/ \ /  \
#    3---2-1--8
nbors = {
    0:[], 1:[2,8], 2:[1,3,4,5], 3:[2,5,9], 4:[2,5,7], 5:[2,3,4,9], 6:[9,],
    7:[4,], 8:[1,4], 9:[3,5,6]
}
print(find_all_rings_dfs(nbors))

# 6--5-8---3-----2
#    |/ \/    \
#    1  /\--4--7
#     \|   /
#      9--/
nbors = {
    0:[], 1:[5,8,9], 2:[3,], 3:[2,7,8,9], 4:[7,8,9], 5:[1,6,8], 6:[5,],
    7:[3,4], 8:[1,3,4,5], 9:[1,3,4]
}
print(find_all_rings_dfs(nbors))



