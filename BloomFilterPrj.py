'''
Bloom Filter
@author: Gabriel A. Granell Jimenez
@date: December 11, 2022
The Goal of project one is to create a Bloom Filter. based on multiple hashes.
'''
import hashlib
import sys
import math

def bitarray(n):
    #Create a bit array of size n
    return [False for _ in range(n)]

#Bloom filter getters:

#Size of the Bloom Filter
def get_size(n, p):
    m = -(n * math.log(p)) / (math.log(2) ** 2)
    return int(m)

#Number of hash functions
def get_hash_count(m, n):
    k = (m / n) * math.log(2)
    return int(k)

# Bloom filter function:
def BloomFilter(n, p):
    #Bloom filter functions
    
    size = get_size(n, p)
    hash_count = get_hash_count(size, n)
    bit_array = bitarray(size)
    for i in range(len(bit_array)):
        bit_array[i] = 0

    # Add an item to the Bloom filter
    def add(item):
        for i in range(hash_count):
            digest = get_hash(item, i)
            bit_array[digest] = True

    # Check if an item is in the Bloom filter
    def check(item):
        for i in range(hash_count):
            digest = get_hash(item, i)
            if bit_array[digest] == False:
                return False
        return True

    #Digest=hash function -> this is to not have to do alot of things to get hash, due to mmh3 not working in Moodle
    def get_hash(item, seed):#Basically digest = hash(item, i) % self.size, it serves as the hash function
        digest = hashlib.sha256()
        # Encode the string to bytes
        digest.update(str(seed).encode() + item.encode())
        return int(digest.hexdigest(), 16) % size

    return add, check
  
def main():
    if len(sys.argv) > 1:
        #Number of emails in the first file
        with open(sys.argv[1]) as file:
            n = sum(1 for line in file)
            
        #Create the Bloom filter
        add, check = BloomFilter(n, 0.0000001)
        # Add the emails from the first file to the Bloom filter
        
        with open(sys.argv[1]) as file:
            next(file)
            for line in file:
                add(line.strip())
                
        # Checking emails from the second file against the created Bloom Filter
        with open(sys.argv[2]) as file:
            next(file)
            for line in file:
                if check(line.strip()):
                    print(line.strip() + ",Probably in the DB")#Is probably in the DB
                else:
                    print(line.strip() + ",Not in the DB")#Is not in the DB

if __name__ == '__main__':
  main()