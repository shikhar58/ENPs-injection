

#ftcs scheme for mass transport  through diffusion and advection
#this is coupled with colloid deposition
#extreme fine x
def mission():
    D=float(e4.get())
    u=float(e3.get())
    ka=float(e8.get())
    kd=float(e9.get())
    rho_b=float(e6.get())
    n=float(e5.get())
    stop=int(e10.get())
    dt=0.01
    dx=0.01
    L=float(e1.get())
    Ltimes=int(L/dx)
    tmax=float(e2.get())
    ttimes=int(tmax/dt)
    c=[float(e7.get())]
    for i in range(1,Ltimes):
        c.append(0)
    c.append(c[Ltimes-1])
    y=[0]
    t=[0]
    s=[0]
    for i in range(1,Ltimes):
        s.append(0)
    s.append(s[Ltimes-1])
    for j in range(1,ttimes):
        if t[j-1]<stop:
            cnew=[0.1]
        else:
            cnew=[0]
        snew=[0]
        for i in range(1,Ltimes):
            p=((D*dt/dx**2))*c[i+1]+(1-u*dt/(n*dx)-2*D*dt/(dx**2)-ka*dt)*c[i]+(D*dt/(dx**2)+u*dt/(n*dx))*c[i-1]+rho_b*kd*dt*s[i]/n
            cnew.append(p)
            q=(n*ka*dt/rho_b)*cnew[i]+(1-kd*dt)*s[i]
            snew.append(q)
        cnew.append(cnew[Ltimes-1])
        snew.append(snew[Ltimes-1])
        y.append(cnew[Ltimes])
        t.append(t[j-1]+dt)
        c=cnew
        s=snew
    import matplotlib.pyplot as plt
    plt.plot(t, y, 'ro')
    #plt.axis([0, 6, 0, 20])
    plt.show()
	
from tkinter import *

master = Tk()
Label(master, text="Length of the column").grid(row=0)
Label(master, text="time period").grid(row=1)
Label(master, text="darcy velocity").grid(row=2)
Label(master, text="dispersivity").grid(row=3)
Label(master, text="porosity").grid(row=4)
Label(master, text="bulk density").grid(row=5)
Label(master, text="Dirchlet boundary concentration").grid(row=6)
Label(master, text="attachment coefficient").grid(row=7)
Label(master, text="detatchment coefficient").grid(row=8)
Label(master, text="time for injection to stop").grid(row=9)


master.minsize(300, 300)
w = Label(master, text="NMN",font=("Helvetica", 32)).grid(row=12, column=1)
e1 = Entry(master)
e2 = Entry(master)
e3 = Entry(master)
e4 = Entry(master)
e5 = Entry(master)
e6 = Entry(master)
e7 = Entry(master)
e8 = Entry(master)
e9 = Entry(master)
e10 = Entry(master)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)
e3.grid(row=2, column=1)
e4.grid(row=3, column=1)
e5.grid(row=4, column=1)
e6.grid(row=5, column=1)
e7.grid(row=6, column=1)
e8.grid(row=7, column=1)
e9.grid(row=8, column=1)
e10.grid(row=9, column=1)


#Button(master, text='Quit', command=master.quit).grid(row=3, column=0, sticky=W, pady=4)
Button(master, text='Run', command=mission).grid(row=10, column=1, sticky=W, pady=4)

mainloop()