import pyarrow.parquet as pq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parquet_file = pq.ParquetFile('/eos/experiment/wcte/wcte_tests/mPMT_periodic_readout/2024-11-24/periodic_readout_20241124210744_2_waveforms.parquet')

batch_size = 10000
batch_id = 0

amplitude = []
digiQ = []
digiT = []
pulses_nbins = 100
pulses_width = 0.5
nPMTs_max = 2000
pulses = [[[] for i in range(pulses_nbins)] for j in range(nPMTs_max)]
dt = 8
fDigiTimeOffset = -30.1

pmt_dict = {}
# key = (1, 2)
# if key in d  #True if exist
# d.get(key)  # Return value if exist, None if not exist

for batch in parquet_file.iter_batches(batch_size=batch_size):
    print("batch ",batch_id)
    batch_id += 1
    v = batch.to_pandas()
    waveforms = np.array(list(v["samples"].values))
    card_id = np.array(list(v["card_id"].values))
    chan = np.array(list(v["chan"].values))
    #res=np.array(list(v["samples"].values)).max(axis=1)
    res_idx=waveforms.argmax(axis=1)
    #amplitude.extend(list(res))
    for i in range(len(res_idx)):
        if waveforms[i][res_idx[i]]>20 and res_idx[i]>=5 and res_idx[i]<=251 and np.sum(waveforms[i][res_idx[i]-2:res_idx[i]+5])>40:
            key = (card_id[i],chan[i])
            pmt_id = -1
            if key in pmt_dict:
                pmt_id = pmt_dict[key]
            else:
                pmt_id = len(pmt_dict)
                pmt_dict[key] = pmt_id
            #amplitude.append(waveforms[i][res_idx[i]])
            #digiQ.append(np.sum(waveforms[i][res_idx[i]-5:res_idx[i]+3]))
            charge = 0
            if waveforms[i][res_idx[i]+2]>0:
                charge = np.sum(waveforms[i][res_idx[i]-5:res_idx[i]+3])
            else:
                charge = np.sum(waveforms[i][res_idx[i]-5:res_idx[i]+2])
            val1 = waveforms[i][res_idx[i]] - waveforms[i][res_idx[i]-1]
            val2 = waveforms[i][res_idx[i]+1] - waveforms[i][res_idx[i]]
            digiT.append(8*val1/(val1-val2))
            time = dt*(res_idx[i]+val1/(val1-val2))+fDigiTimeOffset
            #print("time = ",time)
            #for j in range(res_idx[i]-3,min(res_idx[i]+18,256)):
            for j in range(res_idx[i]-3,res_idx[i]+3):
                idx = int((j*dt-time)/pulses_width)
                if idx>=0 and idx<pulses_nbins:
                    val = waveforms[i][j]/charge
                    if val>-1 and val<1:
                        pulses[pmt_id][idx].append(val)
                        # if idx>0:
                        #     pulses[pmt_id][idx-1].append(val)
                        # if idx+1<pulses_nbins:
                        #     pulses[pmt_id][idx+1].append(val)
    # if batch_id > 10:
    #     break
    # data_dict = batch.to_pydict()
    # for i in range(len(data_dict['samples'])):
    #     amplitude.append(max(data_dict['samples'][i]))

#n, bins, patches = plt.hist(amplitude, bins=range(6,300))
#n, bins, patches = plt.hist(digiQ, bins=range(0,600))
n, bins, patches = plt.hist(digiT,bins=16,range=[0,8])
print(n)
print(bins)
plt.savefig("fig/digiT.pdf")
plt.clf()
#plt.show()

for (card_id, chan), value in pmt_dict.items():
    print(f"fig/{card_id}_{chan}.pdf")
    x=[]
    y=[]
    yerr=[]
    yerr_frac=[]
    for i in range(pulses_nbins):
        if len(pulses[value][i])>0:
            x.append(i*pulses_width)
            y.append(np.mean(pulses[value][i]))
            yerr.append(np.std(pulses[value][i]))
            yerr_frac.append(abs(np.std(pulses[value][i])/np.mean(pulses[value][i])))
    plt.errorbar(x,y,yerr=yerr)
    plt.plot(x,yerr_frac,'o-')
    plt.ylim(-0.15, 0.5)
    plt.grid(visible=True,which='both')
    plt.savefig(f"fig/{card_id}_{chan}.pdf")
    plt.clf()
