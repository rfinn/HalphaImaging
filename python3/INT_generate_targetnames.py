#!/usr/bin/python
"""
Quick script for generating list of target names to put into theli

"""

startN = 16
endN = 33
s = ''
for i in range(startN,endN+1):
    s += f'target-r_{i} '
print(s)

print()
print()


startN = 2
endN = 33
s = ''
for i in range(startN,endN+1):
    s += f'target-Halpha_{i} '
print(s)
