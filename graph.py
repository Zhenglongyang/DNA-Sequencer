
class StringHasher:

    def __init__(self, hash_factor, compression_factor, compression_summand, buckets):
        self.hash_factor = hash_factor
        self.compression_factor = compression_factor
        self.compression_summand = compression_summand
        self.buckets = buckets

    def hash(self, string):
        raw_hash = polynomial_accumulation_hash(self.hash_factor, string, ord)
        return compress(self.compression_factor, self.compression_summand, self.buckets, raw_hash)

def polynomial_accumulation_hash(factor, lst, transform):
    total = 0
    length = len(lst)
    for i in range(length):
        total *= factor
        char_value = transform(lst[length - i - 1])
        total += char_value
    return total

def compress(a, b, N, y):
    return (a * y + b) % N

#Hashtable avec utilisation de polynomial accumulation hash

class HashTable:
    HASH_FACTOR = 41
    size = 0

    def __init__(self, size):
        size = next_prime(size * 2)
        compression_factor = size - 1
        compression_summand = compression_factor - 1
        self.hasher = StringHasher(self.HASH_FACTOR, compression_factor, compression_summand, size)
        self.buckets = []

        for i in range(size):
            self.buckets.append([])

    def get_size(self):
        return self.size

    def get_capacity(self):
        return len(self.buckets)

    def __getitem__(self, key):
        hash = self.hasher.hash(key)

        kv_pairs = self.buckets[hash]
        for kvp in kv_pairs:
            if kvp[0] == key:
                return kvp[1]
        return None

    def __setitem__(self, key, value):
        hash = self.hasher.hash(key)

        kv_pairs = self.buckets[hash]
        for i in range(len(kv_pairs)):
            kvp = kv_pairs[i]
            if kvp[0] == key:
                # Item already exists, so we replace the value
                kv_pairs[i] = (key, value)
                return
        # Item doesn't exist, so insert it in the bucket and increment size
        kv_pairs.append((key, value))
        self.size += 1

    def __delitem__(self, key):
        hash = self.hasher.hash(key)

        kv_pairs = self.buckets[hash]
        for i in range(len(kv_pairs)):
            kvp = kv_pairs[i]
            if kvp[0] == key:
                del kv_pairs[i]
                self.size -= 1

    def __contains__(self, item):
        return not self[item] is None

    def __iter__(self):
        for bucket in self.buckets:
            for kvp in bucket:
                yield kvp[1]

    def resize(self, new_buckets):
        new_table = HashTable(new_buckets)

        for bucket in self.buckets:
            for kvp in bucket:
                new_table[kvp[0]] = kvp[1]
        self = new_table

    def __str__(self):
        string = "["
        for bucket in self.buckets:
            for kvp in bucket:
                key = kvp[0]
                value = kvp[1]
                string += "(" + str(key) + ", " + str(value) + "), "
        string += "]"
        return string


def is_prime(n):
    """
    Checks division by factors up to sqrt(n) for primality.
    """
    factor = 2
    while factor * factor <= n:
        mod = n % factor
        if mod == 0:
            return False
        factor += 1
    return True


def next_prime(n):
    """
    Returns the least prime number equal to or greater than n.
    """
    while not is_prime(n):
        n += 1
    return n

class DeBrujinGraph:

    def __init__(self, nodes, k=3):
        self.alphabet = ["A","T","C","G"]
        #self.table = HashTable(len(nodes))
        self.table = HashTable(len(nodes))
        for i in range(len(nodes)):
           self.table.__setitem__(nodes[i],nodes[i])
        # initialise la structure de donnees

        #for kmer in generate(nodes,k):
        #    self.table.__setitem__(kmer, kmer)

    # determine si le graphe de Brujin contient le noeud N
    def __contains__(self, N):
        return self.table.__contains__(N)

    # retourne un iterable sur les noeuds du graphe
    def __iter__(self):
        return self.nodes()


    # retourne le facteur de charge de la table de hachage sous-jacente
    def load_factor(self) -> object:
        return self.table.get_size()/self.table.get_capacity()

    # ajoute le noeud N au graphe
    def add(self, N):
        self.table.hasher.hash(N)
        # ajoute le noeud N au graphe

    # enleve le noeud N du graphe
    def remove(self, N):
        self.table.__delitem__(N)
        # enleve le noeud N du graphe

    # retourne un iterable sur les noeuds du graphe
    def nodes(self):
        return self.table.__iter__()

    # retourne tous les predecesseur du noeud N
    def predecessors(self, N):
        for i in self.alphabet:
            candidate = i + N[:-1]
            if self.__contains__(candidate):
                yield candidate


    # retourne tous les successeurs du noeud N'
    def successors(self, N):
        for i in self.alphabet:
            candidate= N[1:] + i
            if self.__contains__(candidate):
                yield candidate

    # trouve le noeud qui nont pas de predecesseurs
    def find_start(self):
        for node in self.nodes():
            if not any(self.predecessors(node)):
                yield node
