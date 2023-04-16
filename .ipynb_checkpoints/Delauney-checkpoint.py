import numpy as np
import matplotlib.pyplot as plt

def solve_kepler(M, e, tol):
    
    """
    Solves the kepler equation:
    M = E - e*sin(E)
    
    - M is mean anomaly (input, radians)
    - E is eccentric anomaly (output, radians)
    - e is the eccentricity
    - g(E) = E - e*sin(E) - M
    - g'(E) = 1 - e*cos(E)
    
    Method used: Newtons
    - E_ip1 = Ei - (g(Ei)/g'(Ei))
    - tol: tolerance to stop solving
    """
    
    # initial guess
    E0 = M
    g = E0 - e*np.sin(E0) - M
    
    # iteration number 
    p = 0
    
    Ei = E0
    
    print(np.abs(g))
    
    while np.abs(g) > tol:
        
        p = p+1
        
        # compute g(E) and g'(E)
        g = Ei - e*np.sin(Ei) - M
        g_der = 1 - e*np.cos(Ei)
        
        # update E
        Eip1 = Ei - (g/g_der)
        
        # Compute g at updated value
        g = Eip1 - e*np.sin(Eip1) - M
        
        Ei = Eip1
        
        
    print(f'solution converged in {p} iterations, E = {Eip1}')
    
    return Eip1


e = 0.5
t = np.linspace(1e-5,100,1000)
E = []

for i in range(len(t)):
    
    M = t[i]
    Et = solve_kepler(M=M, e=e, tol=1e-7)

    E.append(Et)
    

# Distance from focus    
r = 1 - e*np.cos(E)

# True anomaly

E = np.array(E)
theta = 2*np.arctan(np.sqrt(3)*np.tan(E/2))

# # Plot r vs t
# plt.figure()
# plt.plot(t, r)
# plt.xlabel('time')
# plt.ylabel('r')
# plt.savefig('Plots/r_vs_t.png')
# plt.show()

x = r*np.cos(theta)
y = r*np.sin(theta)


plt.figure()
plt.plot(0,0,'.')
plt.plot([min(x), max(x)], [0,0])
plt.plot(x,y,'-.')
plt.savefig('Plots/ellipse.png')
plt.show()


# Plot RAAN vs t
w_0 = 0.01
W_0 = 30

W = w_0*t + W_0

plt.figure()
plt.plot(t, W)
plt.xlabel('time')
plt.ylabel('RAAN')
plt.savefig('Plots/RAAN.png')
plt.show()