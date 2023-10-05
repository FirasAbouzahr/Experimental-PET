from PETheader import *

# for obtaining energy/timing resolutions
def gaussian(x,A,mu,sig):
    return A * np.exp(-((x-mu)/sig)**2)

def SingleChannelEnergyResponse(df,channelID,bins,sigma_cut=2.5):
    if channelID in np.unique(df.ChannelIDL):
        df = df[df.ChannelIDL == channelID]
        data = df.ChargeL
    else:
        df = df[df.ChannelIDR == channelID]
        data = df.ChargeR
    
    fig,ax = plt.subplots(figsize = (10,7))
    y,x,_ = plt.hist(data,bins = bins)
    bincenters = np.array([0.5 * (x[i] + x[i+1]) for i in range(len(x) - 1)])
    
    A = np.max(y)
    mu = x[np.where(y == A)[0][0]]
    guess = [A,mu,1]
    
    try:
        p,c = curve_fit(gaussian,bincenters,y,p0=guess)
        xspace = np.linspace(p[1]-2.5*p[2],p[1]+2.5*p[2])
        ax.plot(xspace,gaussian(xspace,*p),color = 'red',linestyle='dashed')
        FWHM = abs(2.355 * p[2])
        photopeakcut = p[1]-sigma_cut*p[2]
        eres = FWHM / p[1] * 100
    except:
        print('Fit-Failed')
        eres = None
        photopeakcut = None
    
    return eres,photopeakcut,p


def getChannelPairs(df,threshold):
    channelData = df.drop(['TimeL','ChargeL','TimeR','ChargeR'],axis=1)
    uniqueChannelPairs,Occurences = np.unique(channelData.to_numpy(),axis = 0,return_counts = True)
    mostActiveChannels = np.where(Occurences >= threshold)[0]
    uniqueChannelPairs = uniqueChannelPairs[mostActiveChannels]
    return uniqueChannelPairs


def getCoincidenceTimeDiffs(df,IDL,IDR,bins,photocut = False,photopeakcuts=[0,0]):
    df_coinc = df[df.ChannelIDL == IDL]
    df_coinc = df_coinc[df_coinc.ChannelIDR == IDR]
    
    # cut data such that it only encompasses data in the photopeak
    if photocut == True:
        df_coinc = df_coinc[df_coinc.ChargeL >= photopeakcuts[0]]
        df_coinc = df_coinc[df_coinc.ChargeR >= photopeakcuts[1]]
    
    timeDiffs = []
    for timeL,timeR in zip(df_coinc.TimeL,df_coinc.TimeR):
        timeDiffs.append(timeL - timeR)
    
    fig,ax = plt.subplots(figsize = (10,7))
    y,x,_ = plt.hist(timeDiffs,bins = bins)
    bincenters = np.array([0.5 * (x[i] + x[i+1]) for i in range(len(x) - 1)])
    
    A = np.max(y)
    mu = x[np.where(y == A)[0][0]]
    guess = [A,mu,np.std(timeDiffs)]
    
    try:
        p,c = curve_fit(gaussian,bincenters,y,p0=guess)
        xspace = np.linspace(p[1]-2.5*p[2],p[1]+2.5*p[2])
        ax.plot(xspace,gaussian(xspace,*p),color = 'red',linestyle='dashed')
        CTR = abs(2.355 * p[2])
    except:
        print('Fit-Failed')
        CTR = None
        
    return CTR,p
    
