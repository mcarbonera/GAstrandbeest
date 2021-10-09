#%% Exemplo de simulação de pêndulo
from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from collections import deque

G = 9.8  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum 1 in m
L2 = 1.0  # length of pendulum 2 in m
L = L1 + L2  # maximal length of the combined pendulum
M1 = 1.0  # mass of pendulum 1 in kg
M2 = 1.0  # mass of pendulum 2 in kg
t_stop = 5  # how many seconds to simulate
history_len = 500  # how many trajectory points to display

def derivs(state, t):

    dydx = np.zeros_like(state)
    dydx[0] = state[1]

    delta = state[2] - state[0]
    den1 = (M1+M2) * L1 - M2 * L1 * cos(delta) * cos(delta)
    dydx[1] = ((M2 * L1 * state[1] * state[1] * sin(delta) * cos(delta)
                + M2 * G * sin(state[2]) * cos(delta)
                + M2 * L2 * state[3] * state[3] * sin(delta)
                - (M1+M2) * G * sin(state[0]))
               / den1)

    dydx[2] = state[3]

    den2 = (L2/L1) * den1
    dydx[3] = ((- M2 * L2 * state[3] * state[3] * sin(delta) * cos(delta)
                + (M1+M2) * G * sin(state[0]) * cos(delta)
                - (M1+M2) * L1 * state[1] * state[1] * sin(delta)
                - (M1+M2) * G * sin(state[2]))
               / den2)

    return dydx

# create a time array from 0..t_stop sampled at 0.02 second steps
dt = 0.02
t = np.arange(0, t_stop, dt)

# th1 and th2 are the initial angles (degrees)
# w10 and w20 are the initial angular velocities (degrees per second)
th1 = 120.0
w1 = 0.0
th2 = -10.0
w2 = 0.0

# initial state
state = np.radians([th1, w1, th2, w2])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, state, t)

x1 = L1*sin(y[:, 0])
y1 = -L1*cos(y[:, 0])

x2 = L2*sin(y[:, 2]) + x1
y2 = -L2*cos(y[:, 2]) + y1

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, 1.))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], ',-', lw=1)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i*dPhi))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(y), interval=dt*1000, blit=True)
plt.show()

#%% Strandbeests
from numpy import sin, cos, sqrt
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from collections import deque

constA = 1
constB = 1
#constC = 1.3
constC = 1.5
constD = 1
constE = 1
constF = 1
#constG = 1.3
constG = 1.5
constH = 1
constI = 1
constJ = 1
constK = 1
constL = 0.5
constM = 0.5

def definirParametros():
    return {"a": constA,
            "b": constB,
            "c": constC,
            "d": constD,
            "e": constE,
            "f": constF,
            "g": constG,
            "h": constH,
            "i": constI,
            "j": constJ,
            "k": constK,
            "l": constL,
            "m": constM}
    
parametros = definirParametros()

def erroDistancia(P1, P2, d):
    dCalculada = sqrt((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2)
    return dCalculada - d

def calculaErros(Pi, Po, Pa, Pb, Pc, Pd):
    erroA = erroDistancia(Pa, Pb, parametros["a"])
    erroB = erroDistancia(Pa, [0,0], parametros["b"])
    erroC = erroDistancia(Pa, Pi, parametros["c"])
    erroD = erroDistancia(Pb, [0,0], parametros["d"])
    erroE = erroDistancia(Pb, Pc, parametros["e"])
    erroF = erroDistancia(Pd, [0,0], parametros["f"])
    erroG = erroDistancia(Pd, Pi, parametros["g"])
    erroH = erroDistancia(Pc, Pd, parametros["h"])
    erroI = erroDistancia(Pc, Po, parametros["i"])
    erroJ = erroDistancia(Pd, Po, parametros["j"])
    
    return [erroA, erroB, erroC, erroD, erroE, erroF, erroG, erroH, erroI, erroJ]

## Derivadas dos Estados
def cinematicaDiretaPi(phi):
    return [parametros["k"] + parametros["m"] * cos(phi),
            parametros["l"] + parametros["m"] * sin(phi)]
    #return [-parametros["k"] - parametros["m"] * cos(phi),
    #        -parametros["l"] - parametros["m"] * sin(phi)]

def cinematicaDiretaPa(Pi):
    tempA = -Pi[0]/Pi[1];
    tempB = (Pi[0]**2 + Pi[1]**2 + parametros["b"]**2 - parametros["c"]**2)/ \
            (2*Pi[1])
    Pax = (-tempA*tempB - sqrt((tempA**2)*(tempB**2) - \
           ((tempA**2)+1)*((tempB**2) - (parametros["b"]**2)))) / \
            ((tempA**2)+1)
    Pay = tempA * Pax + tempB    
    return [Pax, Pay]

def cinematicaDiretaPd(Pi):
    tempC = -Pi[0]/Pi[1];
    tempD = (Pi[0]**2 + Pi[1]**2 + parametros["f"]**2 - parametros["g"]**2)/ \
            (2*Pi[1])
    Pdx = (-tempC*tempD + sqrt((tempC**2)*(tempD**2) - \
           ((tempC**2)+1)*((tempD**2) - (parametros["f"]**2)))) / \
            ((tempC**2)+1)
    Pdy = tempC * Pdx + tempD
    return [Pdx, Pdy]
    
def cinematicaDiretaPb(Pa):
    tempE = -Pa[0]/Pa[1];
    tempF = (parametros["d"]**2 + parametros["b"]**2 - parametros["a"]**2)/ \
            (2*Pa[1])
    Pbx = (-tempE*tempF - sqrt((tempE**2)*(tempF**2) - \
           ((tempE**2)+1)*((tempF**2) - (parametros["d"]**2)))) / \
            ((tempE**2)+1)
    Pby = tempE * Pbx + tempF
    return [Pbx, Pby]

def Bhaskara(A, B, C, plusminus):
    return (-B + (plusminus)*sqrt(B**2 - 4*A*C)) / (2*A)
    
def cinematicaDiretaPc(Pb, Pd):    
    tempParametros = (parametros["e"]**2 - parametros["d"]**2 + \
        parametros["f"]**2 - parametros["h"]**2)
    tempPcxParcial1 = ((Pb[0] - Pd[0])**2 + (Pb[1] - Pd[1])**2)
    tempPcxParcial2 = (Pb[0] - Pd[0])*(2*Pb[1]*(Pb[1] - Pd[1]) + \
        tempParametros) - 2*Pb[0]*((Pb[1] - Pd[1])**2)
    tempPcxParcial3 = Pb[1]*(Pb[1] - Pd[1])*(tempParametros) - \
        ((Pb[1] - Pd[1])**2)*(parametros["e"]**2 - parametros["d"]**2) + \
            ((tempParametros)/(2))**2
    ctteSolucao = -1
    Pcx = Bhaskara(tempPcxParcial1, tempPcxParcial2, tempPcxParcial3, ctteSolucao)
    Pcy = Pb[1] - sqrt(Pb[1]**2 - Pcx**2 + 2*Pb[0]*Pcx - parametros["d"]**2 +\
                       parametros["e"]**2)
    return [Pcx, Pcy]
    
def calculaVetorErro(PoxVet, Poy, Pc, Pd, dc, dd):
    Pox1 = PoxVet[0]
    Pox2 = PoxVet[1]    
    return [
    sqrt(erroDistancia([Pox1, Poy], Pc, dc)**2 + \
         erroDistancia([Pox1, Poy], Pd, dd)**2),
    sqrt(erroDistancia([Pox2, Poy], Pc, dc)**2 + \
         erroDistancia([Pox2, Poy], Pd, dd)**2)
    ]
    
def cinematicaDiretaPo(Pb, Pc, Pd):
    varTemp = Pc[0]**2 - Pc[1]**2 - 2*Pc[0]*Pd[0] + parametros["i"]**2 + \
        parametros["f"]**2 - parametros["j"]**2
    tempPoxParcial1 = (Pc[1]-Pd[1])**2 + (Pd[0]-Pc[0])**2
    tempPoxParcial2 = (Pc[1] - Pd[1])*varTemp - 2*Pc[1]*((Pd[0]-Pc[0])**2)
    tempPoxParcial3 = (Pc[1]**2 - parametros["i"]**2)*((Pd[0]-Pc[0])**2) + \
        (varTemp/2)**2
    ctteSolucao = -1
    Poy = Bhaskara(tempPoxParcial1, tempPoxParcial2, tempPoxParcial3, ctteSolucao)
    
    cttePox = np.array([-1, 1])
    PoxVet = Pc[0] + cttePox * \
        sqrt(parametros["i"]**2 + 2*Pc[1]*Poy - Pc[1]**2 - Poy**2)
    
    Pox = 0
    erroVet = calculaVetorErro(PoxVet, Poy, Pc, Pd, parametros["i"], 
                               parametros["j"])
    if(erroVet[0] < erroVet[1]):
        Pox = PoxVet[0]
    else:
        Pox = PoxVet[1]
    
    #if(np.abs(erroDistancia([Pox, Poy], Pd, parametros["j"])) > 0.001):
    #    print('Pobrema: ', erroDistancia([Pox, Poy], Pd, parametros["j"]),
    #          varTemp, tempPoxParcial1, tempPoxParcial2, 
    #          tempPoxParcial3, Pox, Poy)
    return [Pox, Poy]


def calculaEstado(phi):
    Pi = cinematicaDiretaPi(phi)
    Pa = cinematicaDiretaPa(Pi)
    Pd = cinematicaDiretaPd(Pi)
    Pb = cinematicaDiretaPb(Pa)
    Pc = cinematicaDiretaPc(Pb, Pd)
    Po = cinematicaDiretaPo(Pb, Pc, Pd)
    return {
        "Pi": Pi,
        "Pa": Pa,
        "Pd": Pd,
        "Pb": Pb,
        "Pc": Pc,
        "Po": Po,
        "raw": Pi + Pa + Pd + Pb + Pc + Po
    }

# INICIAR SIMULACAO.

# Variável phi:
phi_stop = 2*np.pi

# phi0 é o angulo inicial (graus)
# omega0 é a velocidade angular. (graus por segundo)
phi0 = 0.0
omega0 = np.radians(90.0)

dPhi = 0.05/np.pi
phi = np.arange(phi0, phi_stop, dPhi)
history_len = len(phi) # how many trajectory points to display

y = []
for phiValue in phi:
    estadoTemp = calculaEstado(phiValue)
    y.append(estadoTemp["raw"])

## Animação
fig = plt.figure(figsize=(5, 4))
LIM = 2
offsetY = -0.5
ax = fig.add_subplot(autoscale_on=False, xlim=(-LIM, LIM), 
                     ylim=(-LIM+offsetY, LIM+offsetY))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], ',-', lw=1)
phi_template = 'phi = %.1f'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)    
    
def animate(i):
    estadoAtual = y[i]
    Pi = [estadoAtual[0], estadoAtual[1]]
    Pa = [estadoAtual[2], estadoAtual[3]]
    Pd = [estadoAtual[4], estadoAtual[5]]
    Pb = [estadoAtual[6], estadoAtual[7]]
    Pc = [estadoAtual[8], estadoAtual[9]]
    Po = [estadoAtual[10],estadoAtual[11]]
    
    thisx = [0, Pi[0], Pa[0], Pd[0], Pb[0], Pc[0], Po[0]]
    thisy = [0, Pi[1], Pa[1], Pd[1], Pb[1], Pc[1], Po[1]]
    
    #print(calculaErros(Pi, Po, Pa, Pb, Pc, Pd))
    
    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[6])
    history_y.appendleft(thisy[6])

    lineX = [0, Pd[0], Pc[0], Pb[0], Pc[0], Po[0], Pd[0], Pi[0], Pa[0], 0, Pb[0], Pa[0]]
    lineY = [0, Pd[1], Pc[1], Pb[1], Pc[1], Po[1], Pd[1], Pi[1], Pa[1], 0, Pb[1], Pa[1]]
    line.set_data(lineX, lineY)
    
    trace.set_data(history_x, history_y)
    time_text.set_text(phi_template % (i*dPhi))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(phi), interval=dPhi*1000, blit=True)
plt.show()