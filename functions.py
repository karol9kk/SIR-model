# download modsim.py if necessary

from os.path import basename, exists

def download(url):
    filename = basename(url)
    if not exists(filename):
        from urllib.request import urlretrieve
        local, _ = urlretrieve(url, filename)
        print('Downloaded ' + local)
    
download('https://raw.githubusercontent.com/AllenDowney/' +
         'ModSimPy/master/modsim.py')

# import from modsim
from modsim import *

def make_system_number(beta, gamma,time_predicted,init,vaccine_rate):
    
    return System(init=init, t_end=time_predicted,
                  beta=beta, gamma=gamma,vaccine_rate=vaccine_rate)

def make_system(beta, gamma,time_predicted,init,vaccine_rate):
    
    init /= init.sum()

    return System(init=init, t_end=time_predicted,
                  beta=beta, gamma=gamma,vaccine_rate=vaccine_rate)


def update_func(t, state, system):
    
    s, i, r,v = state.s, state.i, state.r,state.v

    
    vaccine_rate=system.vaccine_rate

    infected = system.beta * s* i    
    recovered = system.gamma * i
    
    vacinated=(s+r)*vaccine_rate
    vacinated_s=s*vaccine_rate
    vacinated_r=r*vaccine_rate
    
    
    recovered_not_immune=r*0.5
    
    s -= (infected+vacinated_s-recovered_not_immune)
    i += (infected - recovered)
    r += (recovered- vacinated_r-recovered_not_immune)
    v += vacinated
    
    return State(s=s, i=i, r=r,v=v)


def run_simulation(system, update_func):
    frame = TimeFrame(columns=system.init.index)
    frame.loc[0] = system.init
    
    for t in range(0, system.t_end):
        frame.loc[t+1] = update_func(t, frame.loc[t], system)
        
    
    return frame



def plot_results(S, I, R,V):
    
    S.plot(style='--', label='Susceptible')
    I.plot(style='-', label='Infected')
    R.plot(style=':', label='Recovered')
    V.plot(style='-.',label='Vacniated')
    decorate(xlabel='Time (days)',
             ylabel='Fraction of population')
    plt.grid(True)
    





def calc_total_infected(results, system):
    s_0 = results.s[0]
    s_end = results.s[system.t_end]
    return s_0 - s_end

def sweep_immunity(fraction_array,beta, gamma,time_predicted,init):
    sweep = SweepSeries()

    for fraction in fraction_array:
        system = make_system(beta, gamma,time_predicted,init)
        results = run_simulation(system, update_func)
        sweep[fraction] = calc_total_infected(results, system)

    return sweep