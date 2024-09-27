from collections import namedtuple
import time

# Define namedtuple
MyTuple = namedtuple('MyTuple', ['field1', 'field2'])

# Create namedtuple and dictionary
nt = MyTuple(1, 2)
d = {'field1': 1, 'field2': 2}

# Access benchmark
start_time = time.time()
for _ in range(1000000):
    x = nt.field1
    #nt.field1 = 3
namedtuple_time = time.time() - start_time

start_time = time.time()
for _ in range(1000000):
    x = d['field1']
    #d['field1'] = 3
dict_time = time.time() - start_time

print(f"Namedtuple access time: {namedtuple_time}")
print(f"Dictionary access time: {dict_time}")
