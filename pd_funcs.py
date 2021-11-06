def get_neighbors(df,current_id,horizon):
    distance=((df.coordx-df.loc[current_id].coordx)**2+(df.coordy-df.loc[current_id].coordy)**2)**0.5
    neighbor_ids = distance[distance<=horizon].index
    neighbor_ids=neighbor_ids.drop(current_id)
    idist = ((df.loc[neighbor_ids].coordx-df.loc[current_id].coordx)**2+(df.loc[neighbor_ids].coordy-df.loc[current_id].coordy)**2)**0.5
    idist.name = 'idist'
    return neighbor_ids,idist

def compute_nlength(df,current_id,neighbor_ids):
    nlength = (((df.loc[current_id].coordx + df.loc[current_id].dispx)-(df.loc[neighbor_ids].coordx + df.loc[neighbor_ids].dispx))**2+((df.loc[current_id].coordy + df.loc[current_id].dispy)-(df.loc[neighbor_ids].coordy + df.loc[neighbor_ids].dispy))**2)
    nlength.name = 'nlength'
    return nlength
	
def compute_stretch(nlength,idist):
    stretch = nlength/idist
    stretch.name = 'stretch'
    return stretch

def compute_fac_volume_cor(idist,horizon,delta):
    r = delta/2
    fac_volume_corr=pandas.Series(index=idist.index,dtype='float32')
    fac_volume_corr[idist[idist<=horizon-r].index]=1.
    fac_volume_corr[idist[idist>horizon-r][idist[idist>horizon-r]<=horizon+r].index]=(horizon+r)-idist.loc[idist[idist>horizon-r][idist[idist>horizon-r]<=horizon+r].index]/2/r
    fac_volume_corr[idist[idist>horizon+r].index]=0.
    fac_volume_corr.name = 'fac_volume_corr'
    return fac_volume_corr
	
def compute_Lambda(df,current_id,neighbor_ids,idist,nlength):
    ykx=(df.loc[current_id].dispx+df.loc[current_id].coordx)
    yky=(df.loc[current_id].dispy+df.loc[current_id].coordy)
    yjx=(df.loc[neighbor_ids].dispx+df.loc[neighbor_ids].coordx)
    yjy=(df.loc[neighbor_ids].dispy+df.loc[neighbor_ids].coordy)
    xkx=(df.loc[current_id].coordx)
    xky=(df.loc[current_id].coordy)
    xjx=(df.loc[neighbor_ids].coordx)
    xjy=(df.loc[neighbor_ids].coordy)
    Lambda=((yjx-ykx)*(xjx-xkx)+(yjy-yky)*(xjy-xky))/idist/nlength
    Lambda.name = 'Lambda'
    return Lambda

def compute_Dilatation(df,current_id,neighbor_ids,stretch,fac_volume_corr,Lambda,horizon,d):
    Dilatation = horizon*d*Lambda*stretch*df.loc[neighbor_ids].volume * fac_volume_corr
    Dilatation.name = 'Dilatation'
    return Dilatation
	
def compute_SED(df,current_id,neighbor_ids,nlength,idist,fac_volume_corr,horizon,b,Dilatation,a):
    SED = horizon*b*((nlength-idist))**2 /idist * df.loc[neighbor_ids].volume*fac_volume_corr
    SED.name = 'SED'
    return SED

def construct_dfn(idist,nlength,stretch,fac_volume_corr,Lambda,Dilatation,SED):
    dfn = pandas.DataFrame()
    dfn[idist.name]=idist
    dfn[nlength.name]=nlength
    dfn[fac_volume_corr.name]=fac_volume_corr
    dfn[Lambda.name]=Lambda
    dfn[Dilatation.name]=Dilatation
    dfn[SED.name]=SED
    return dfn
	
def set_surface_correction_factors(df,current_id,condition,disp_grad,SED):
    if condition=='uniaxial stretch x':
        df.loc[current_id].D1 = disp_grad / df.loc[current_id].Dilatation
        df.loc[current_id].Dilatation = df.loc[current_id].Dilatation *df.loc[current_id].D1 
    elif condition=='uniaxial stretch y':
        df.loc[current_id].D2 = disp_grad / df.loc[current_id].Dilatation
        df.loc[current_id].Dilatation = df.loc[current_id].Dilatation *df.loc[current_id].D2 
    elif condition=='simple shear in x-y':
        df.loc[current_id].S1 = 0.5*disp_grad**2*mu / SED.sum()
        df.loc[current_id].S2 = 0.5*disp_grad**2*mu / SED.sum()
        df.loc[current_id].SED = SED.sum()*df.loc[current_id].S1

def preprocess(df,condition_list,disp_grad):
    for condition in condition_list:
        set_loading_condition(condition,df,disp_grad)
        for current_id in df.index:
            neighbor_ids,idist=get_neighbors(df,current_id,horizon)
            nlength=compute_nlength(df,current_id,neighbor_ids)
            stretch=compute_stretch(nlength,idist)
            fac_volume_corr=compute_fac_volume_cor(idist,horizon,delta)
            Lambda=compute_Lambda(df,current_id,neighbor_ids,idist,nlength)
            Dilatation=compute_Dilatation(df,current_id,neighbor_ids,stretch,fac_volume_corr,Lambda,horizon,d)
            df.loc[current_id].Dilatation= Dilatation.sum()
            SED=compute_SED(df,current_id,neighbor_ids,nlength,idist,fac_volume_corr,horizon,b,Dilatation,a)
            df.loc[current_id].SED = a*Dilatation.sum()**2 + SED.sum()
            dfn=construct_dfn(idist,nlength,stretch,fac_volume_corr,Lambda,Dilatation,SED)
            set_surface_correction_factors(df,current_id,condition,disp_grad,SED)

def surface_correction_vector(df,current_id,neighbor_ids,idist):
    gdx=(df.loc[current_id].D1 + df.loc[neighbor_ids].D1)*0.5
    gdy=(df.loc[current_id].D2 + df.loc[neighbor_ids].D2)*0.5
    gbx=(df.loc[current_id].S1 + df.loc[neighbor_ids].S1)*0.5
    gby=gbx
    nx=(df.loc[neighbor_ids].coordx-df.loc[current_id].coordx)/idist
    ny=(df.loc[neighbor_ids].coordy-df.loc[current_id].coordy)/idist
    Gd=((nx/gdx)**2+(ny/gdy)**2)**-0.5
    Gb=((nx/gbx)**2+(ny/gby)**2)**-0.5
    return Gd,Gd

def preprocess_with_SCF(df):
    for current_id in df.index:
        neighbor_ids,idist=get_neighbors(df,current_id,horizon)
        Gd, Gb =surface_correction_vector(df,current_id,neighbor_ids,idist)
        nlength=compute_nlength(df,current_id,neighbor_ids)
        stretch=compute_stretch(nlength,idist)
        fac_volume_corr=compute_fac_volume_cor(idist,horizon,delta)
        Lambda=compute_Lambda(df,current_id,neighbor_ids,idist,nlength)
        Dilatation=compute_Dilatation(df,current_id,neighbor_ids,stretch,fac_volume_corr,Lambda,horizon,d)
        Dilatation = Dilatation * Gd
        df.loc[current_id].Dilatation= Dilatation.sum()
        SED=compute_SED(df,current_id,neighbor_ids,nlength,idist,fac_volume_corr,horizon,b,Dilatation,a)
        SED = SED*Gb
        df.loc[current_id].SED = a*Dilatation.sum()**2 + SED.sum()
        dfn=construct_dfn(idist,nlength,stretch,fac_volume_corr,Lambda,Dilatation,SED)

def compute_PD_forces(df,current_id,neighbor_ids,idist,nlength,stretch,fac_volume_corr,Lambda):
    A = 2*horizon*(d*Lambda/idist*(a*df.loc[current_id].Dilatation)+b*(stretch))
    A.name = 'A'
    B = 2*horizon*(d*Lambda/idist*(a*df.loc[neighbor_ids].Dilatation)+b*(stretch))
    delyx =(df.loc[neighbor_ids].dispx +df.loc[neighbor_ids].coordx - (df.loc[current_id].dispx+df.loc[current_id].coordx))
    delyy = (df.loc[neighbor_ids].dispy +df.loc[neighbor_ids].coordy - (df.loc[current_id].dispy+df.loc[current_id].coordy))
    tkjx = 0.5 * A * delyx/nlength
    tkjy = 0.5 * A * delyy/nlength
    tjkx = -0.5 * B * delyx/nlength
    tjky = -0.5 * B * delyy/nlength
    pforcex = ((tkjx-tjkx)*fac_volume_corr*df.loc[neighbor_ids].volume).sum()
    pforcey = ((tkjy-tjky)*fac_volume_corr*df.loc[neighbor_ids].volume).sum()
    return pforcex, pforcey
	
def compute_stable_mass_vector(df,horizon,current_id,neighbor_ids,idist,fac_volume_corr,dt,a,b,safety_factor):
    Kijx = numpy.abs(df.loc[neighbor_ids].coordx-df.loc[current_id].coordx)/idist/idist * 4 * horizon * (0.5 * a * d**2 * horizon / idist * (df.loc[current_id].volume + df.loc[neighbor_ids].volume*fac_volume_corr)+b)
    Kijy = numpy.abs(df.loc[neighbor_ids].coordy-df.loc[current_id].coordy)/idist/idist * 4 * horizon * (0.5 * a * d**2 * horizon / idist * (df.loc[current_id].volume + df.loc[neighbor_ids].volume*fac_volume_corr)+b)
    massvecx = safety_factor*0.25*dt**2*Kijx.sum()
    massvecy = safety_factor*0.25*dt**2*Kijy.sum()
    return massvecx, massvecy

def compute_damping_coeff(df,dt):
    cn=0.
    cn_num=0.
    cn_denom=0.
    for current_id in df.index:
        if (df.loc[current_id].velhalfoldx!=0.):
            Kiix = -(df.loc[current_id].pforcex - df.loc[current_id].pforceoldx)/df.loc[current_id].massvecx/dt/df.loc[current_id].velhalfoldx
            cn_num=cn_num + df.loc[current_id].dispx*Kiix*df.loc[current_id].dispx
        if (df.loc[current_id].velhalfoldy!=0.):
            Kiiy = -(df.loc[current_id].pforcey - df.loc[current_id].pforceoldy)/df.loc[current_id].massvecy/dt/df.loc[current_id].velhalfoldy
            cn_num=cn_num + df.loc[current_id].dispy*Kiiy*df.loc[current_id].dispy
        cn_denom=cn_denom+df.loc[current_id].dispx**2
        cn_denom=cn_denom+df.loc[current_id].dispy**2
    if (cn_denom!=0.):
        if ((cn_num/cn_denom)>0.):
            cn = 2 * (cn_num/cn_denom)**0.5
        else:
            cn=0.
    else:
        cn=0.
    if (cn>2.):
        cn=1.9
    return cn

def apply_ADR(df,tt,dt,cn):
    for current_id in df.index:
        if(tt==1):
            df.loc[current_id].velhalfx = dt * (df.loc[current_id].pforcex+df.loc[current_id].bforcex)/df.loc[current_id].massvecx/2
            df.loc[current_id].velhalfy = dt * (df.loc[current_id].pforcey+df.loc[current_id].bforcey)/df.loc[current_id].massvecy/2
        else:
            df.loc[current_id].velhalfx = ((2.0 - cn * dt)*df.loc[current_id].velhalfoldx+2*dt*(df.loc[current_id].pforcex+df.loc[current_id].bforcex)/df.loc[current_id].massvecx)/(2+cn*dt)
            df.loc[current_id].velhalfy = ((2.0 - cn * dt)*df.loc[current_id].velhalfoldy+2*dt*(df.loc[current_id].pforcey+df.loc[current_id].bforcey)/df.loc[current_id].massvecy)/(2+cn*dt)
        df.loc[current_id].velx = (df.loc[current_id].velhalfoldx+df.loc[current_id].velhalfx)*0.5
        df.loc[current_id].vely = (df.loc[current_id].velhalfoldy+df.loc[current_id].velhalfy)*0.5
        df.loc[current_id].dispx=(df.loc[current_id].dispx+df.loc[current_id].velx*dt)
        df.loc[current_id].dispy=(df.loc[current_id].dispy+df.loc[current_id].vely*dt)
        df.loc[current_id].velhalfoldx = df.loc[current_id].velhalfx
        df.loc[current_id].velhalfoldy = df.loc[current_id].velhalfy
        df.loc[current_id].pforceoldx = df.loc[current_id].pforcex
        df.loc[current_id].pforceoldy = df.loc[current_id].pforcey
	
def iterate(df,ibcs,dt,max_iter,applied):
    set_loading_condition(ibcs,df,applied)
    for tt in range(1,max_iter+dt,dt):
        if numpy.mod(tt,50)==0:
            print(df)
        preprocess_with_SCF(df)
        print(tt)
        for current_id in df.index:
            neighbor_ids,idist=get_neighbors(df,current_id,horizon)
            nlength=compute_nlength(df,current_id,neighbor_ids)
            stretch=compute_stretch(nlength,idist)
            fac_volume_corr=compute_fac_volume_cor(idist,horizon,delta)
            Lambda=compute_Lambda(df,current_id,neighbor_ids,idist,nlength)
            pforcex, pforcey = compute_PD_forces(df,current_id,neighbor_ids,idist,nlength,stretch,fac_volume_corr,Lambda)
            df.loc[current_id].pforcex = pforcex
            df.loc[current_id].pforcey = pforcey
            massvecx, massvecy = compute_stable_mass_vector(df,horizon,current_id,neighbor_ids,idist,fac_volume_corr,dt,a,b,safety_factor=5)
            df.loc[current_id].massvecx = massvecx
            df.loc[current_id].massvecy = massvecy
        cn = compute_damping_coeff(df,dt)
        apply_ADR(df,tt,dt,cn)
