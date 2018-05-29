import networkx as nx
import dwave_networkx as dnx
import minorminer as mm
import sys
from math import log
from random import choice, shuffle


def triangle_g(n, g=None):
    if g is None:
        g = nx.Graph()
    for i in range(n):
        g.add_edge((i, 0), (i, 1))
        for j in range(i):
            g.add_edge((j, 0), (i, 1))
            g.add_edge((j, 1), (i, 1))
    return g


def beta_pass(miner, var_order, beta_schedule, **args):
    args['overlap_bound'] = 9999
    for v in var_order:
        emb = miner.quickpass(
            [v], round_beta=beta_schedule[v], clear_first=False, **args)
        yield emb


def bfs_order(g):
    start = choice(list(g.nodes()))
    order = [start]
    for node, succs in nx.bfs_successors(g, start):
        shuffle(succs)
        order.extend(succs)
    return order


def BetaScheduleHeuristicEmbedding(s, t, beta_schedule):
    miner = mm.miner(list(s.edges()), list(t.edges()))
    # initialization
    for emb in beta_pass(miner, bfs_order(s), beta_schedule):
        pass
    best_emb = emb
    best_key = miner.quality_key(emb)

    for i in range(10):
        for emb in beta_pass(miner, bfs_order(s), beta_schedule):
            key = miner.quality_key(emb)
            if key < best_key:
                best_emb = emb
                best_key = key

    if best_key[1]:
        return {}
#        print "overlaps:", best_key[1]
    else:
        #        print "done"
        return best_emb


g = triangle_g(10)
c = dnx.chimera_graph(8, 8, 2)


for beta_0 in [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]:
    e = g.number_of_edges()+0.
    for beta_add in [lambda d: 0, lambda d: (d**.5)/2, lambda d: log(d)/2., lambda d: d/e]:
        schedule = {v: beta_0 + beta_add(d) for v, d in g.degree()}
        embs = 0
        effort = 0
        scl = 0.
        for _ in range(1000):
            emb = BetaScheduleHeuristicEmbedding(g, c, schedule)
            if emb:
                embs += 1
                scl += sum(len(c) for c in emb.values())
        if embs:
            print (embs, scl/embs),
        else:
            print None,

    print
