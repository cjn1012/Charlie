# Lab 5 Problem 4
def catalan(n):
   if n == 0:
       return 1
   else:
       return ((4*n-2)*catalan(n-1))/(n+1)   

#catalan(100) comes out to be 896519947090131496687170070074100632420837521538745909320

print(catalan(100))

def gcd(m, n):
   if n == 0:
       return m
   else:
       return gcd(n, m%n)

gcd192_108 = gcd(192,108)
print(gcd192_108)

