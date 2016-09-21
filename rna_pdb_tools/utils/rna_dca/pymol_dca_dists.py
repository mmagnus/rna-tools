def get_dists(interactions): # l=([1,2], [3,4])
    for i in interactions:
        a = "////" + str(i[0]) + "/C2'"
        b = "////" + str(i[1]) + "/C2'"
        cmd.distance('d'+str(i[0])+'-'+str(i[1]),"(" + a + ")", "(" + b + ")")# mode, 4)

cmd.extend('get_dists', get_dists)
