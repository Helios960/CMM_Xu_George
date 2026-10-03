import numpy as np  
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--dimpled", action="store_true", help="Is the ball dimpled?")










def Cd(mag):
    if 
    

def main():
    
    
    
    x, y = np.empty(1000), np.empty(1000)
    x[0], y[0] = 0, 0
    vi = 
    thi = 
    vx, vy= np.empty(1000), np.empty(1000)
    vx[0], vy[0] = vi * np.cos(thi), vi * np.sin(thi) 
    time = np.empty(1000)
    mag = np.empty(1000)
    mag[0] = vi
    
    timedelta = 0.01
    rho = 1.225
    A = 
    m = 
    g = 9.81

    i = 0
    
    while y[i] >= 0:
        time[i+1] = (i+1) * timedelta
        x[i+1] = vx[i] * timedelta
        y[i+1] = vy[i] * timedelta
        
        
        mag[i+1] = np.sqrt(vx**2 + vy**2)
        vx[i+1] = (Cd(mag) * rho * A * mag * vx / m) * timedelta
        
        
        vy -= (g + (C * rho * A * mag * vy / m)) * timedelta
        time += timedelta
        i += 1
    