import numpy as np

dxi = (np.pi)/4
theta1 = 1-(1/6)*(dxi**2)
theta2 = -(1/3)*dxi
i_list = np.arange(1,6)


def theta1iplus1(theta1i, theta2i, dxi):
    theta1 = theta1i + theta2i*dxi
    return theta1

def theta2iplus1(i, theta1i, theta2i, dxi):
    theta2 = (((i*dxi)/((i+1)*dxi))**2)*(theta2i - (theta1i**1.5)*dxi)
    return theta2




t1i =[]
t2i = []
xi = []
t1i.append(theta1)
t2i.append(theta2)
xi.append(dxi)

for i in i_list:
    y = theta1
    theta1 = theta1iplus1(theta1, theta2, dxi)
    theta2 = theta2iplus1(i, y, theta2, dxi)
    t1i.append(theta1)
    t2i.append(theta2)
    xi.append((i+1)*dxi)
print(xi)
print(t1i)
print(t2i)

    
