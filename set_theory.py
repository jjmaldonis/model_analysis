""" Set theory operations. """

def setsize(s):
    if not s: return 0.0
    return s[1]-s[0]

def union(*args):
    a = float('inf')
    b = -float('inf')
    for s in args:
        if s[0] < a: a = s[0]
        if s[1] > b: b = s[1]
    return (a, b)

def intersect(*args):
    a = -float('inf')
    b = float('inf')
    for s in args:
        if s[0] > a: a = s[0]
        if s[1] < b: b = s[1]
    if a > b:
        return ()
    return (a, b)

#def complement(s):

def overlap(s1, s2):
    avg = (setsize(s1)+setsize(s2))/2
    i = intersect(s1, s2)
    overlap = setsize(i)/avg*100
    return overlap

def main():
    from math import sqrt
    s1 = (0, 1.3)
    s2 = (0.97, 3.9)

    ico = [555, 504, 594, 546, 625, 547, 592, 551, 540, 595]
    value1 = ico[0]
    value2 = ico[2]
    s1 = (value1-sqrt(value1), value1+sqrt(value1))
    s2 = (value2-sqrt(value2), value2+sqrt(value2))

    print(s1)
    print(s2)
    print(overlap(s1,s2))


if __name__ == '__main__':
    main()
