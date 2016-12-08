import before
import after

flatten = lambda l: [item for sublist in l for item in sublist]

print(len(flatten(flatten(before.test1))))
print(len(flatten(flatten(after.test2))))
