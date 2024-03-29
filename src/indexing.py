import bisect
class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()
        
    def query(self, p):
        kmer = p[:self.k]
        # kmer is the object we are going to search for
        # This assure that we always get the first occurrence of that.
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            # index equals to kmer, then append the location of hte index
            hits.append(self.index[i][1])
            i += 1
        return hits
        
def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        # Verification
        if p[k:] == t[i+k:i+len(p)]:
            offsets.append(i)
    return offsets           
 
t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'
index = Index(t, 2)

test = index.query(p)
print(test)
print(index.query(p))

print(queryIndex(p, t, index))
