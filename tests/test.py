from cycler import cycler
from matplotlib import pyplot
from bioenergetics.util import daphnia

year = '2015'
month = 'June'
site = 'Fall Creek'
total_daphnia = 5020.65
daphnia_size = 1.26
site_data = Site_Data(year, site, month, 0.72, 25, 0.1)
starting_mass = 3
daph_data = Daph_Data(5020.65, 1.26, year, site, month)
max_temp = 10000
min_temp = -1
cust_temp = '{0}_T_{1}_{2}.csv'.format(site, month, year)
elev = 691
pop_site = site

x,y = daphnia.select_rows('data/daphnia-vd.csv', site, month, year)
daph_line, daph_auc = daphnia.compute_curves(x, y, total_daphnia)

def make_plots(self):
    w = 0.5
    temp = 15
    ws = [np.round(x,2) for x in np.arange(0.1,1.0,0.1)]
    ps = np.arange(0,1.0,0.01)
    gs = {}
    temps = [np.round(x,2) for x in np.arange(5,25,2)]
    for w in ws:
        gs[w] = []
        for P in ps:
            (cons,eg,ex,res,sda) = \
                self.compute_bioenergetics(w, temp, P, self.prey,
                                           self.digestibility)
            growth = self.compute_growth(cons, self.prey, self.preyenergy,
                                         eg, ex, sda,res,
                                         self.predatorenergy, w)
            gs[w].append(growth)
    fig,ax = pyplot.subplots()
    colors = [pyplot.get_cmap('inferno')(1.0 * i/len(ws))
              for i in range(len(ws))]
    ax.set_prop_cycle(cycler('color',colors))
    for w,growths in gs.items():
        ax.plot(ps,growths,label=str(w))
    ax.legend()
    pyplot.xlabel('P')
    pyplot.ylabel('growth')
    pyplot.show()


Batch.make_plots = make_plots
FRESH_BATCH = Batch(site_data, starting_mass, daph_data,
                    max_temp, min_temp, cust_temp, elev, pop_site, True)

#FRESH_BATCH.make_plots()

BASE_RESULTS, DAPHNIA_CONSUMED, CONDITION, CONDITION1, DAY_TEMP, NIGHT_TEMP,\
    POPULATION_ESTIMATE = FRESH_BATCH.Run_Batch()
keys = ['StartingMass','consumption','egestion','excretion',
        'P', 'day_depth','night_depth','temps']
rows = np.ceil(len(keys)/2)
for idx,k in enumerate(keys):
    pyplot.subplot(rows,2,idx+1)
    pyplot.plot(BASE_RESULTS[k])
    pyplot.title(k)

depths = np.arange(FRESH_BATCH.depth_min, FRESH_BATCH.depth_max)
fig, ax1 = pyplot.subplots()
ax1.plot(FRESH_BATCH.temp_from_depth(depths), depths, 'orange')
ax2 = ax1.twiny()
ax2.plot(FRESH_BATCH.daphline(depths), depths, 'green')
pyplot.show()
