library(tidyverse)
library(sf)
library(INLA)

rm(list=ls())

#lac <- as_tibble(read.csv("../data/LA_County_ZIP_Codes.csv")) %>%
#    mutate(
#        zip = as.character(ZIPCODE)
#    ) %>%
#    pull(zip)

###########
### Fix ###
###########

## Fix this
dat.org <- as_tibble(read.csv("../data/TEST_20241029_TotalERDaily1618pois.csv")) %>%
    filter(
#        zip %in% lac,
        zip != "90704", ## Catalina
        zip != "95389", ## Yosemite
        zip != "95389", ## West Yosemite
        zip != "95318"  ## El Portal
    ) %>%
    mutate(
      #  cz = sample(1:20,n(),rep=TRUE), ### Fix this
        er = D2Dx1 ### Fix this
    )

## comment out this !!!!! FIX!!!!!!!!!!!!!!!!!!!!
#set.seed(1)
#n12 <- length(dat.org$er[dat.org$er == "<12"])
#dat.org$er[dat.org$er == "<12"] <- "1"#sample(as.character(0:5),n12,rep=TRUE)

############
### Data ###
############

dat <- dat.org %>%
    transmute(
        er = as.integer(er),
        cz = as.integer(FinalClust), #DGG updated 10/31/24
        doy = yday(ymd(date)),
        dow = wday(ymd(date)),
        week = week(ymd(date)),
        year = year(ymd(date)),
        zip = as.character(zip),
        n = as.integer(Population),
        heat = Max_HI_Value
    ) %>%
    drop_na() %>%
    filter(
        heat != 0,
        n != 0L
    ) %>%
    arrange(zip)

cali.shp <- st_read("../data/shapefiles/ca/CA_ZCTA_2020_Final.shp") %>%
    mutate(
        zip = as.character(ZCTA5CE20)
    ) %>%
    filter(
        zip %in% dat$zip
    )

source("../support/make.adj.R")
spat.stuff <- make.adj(cali.shp,dat,"zip","my.zip.adj.txt")
H <- inla.read.graph("my.zip.adj.txt")
cali.shp <- spat.stuff$shapefile

dat.spatial.org <- as_tibble(spat.stuff$datafile) %>%
    
    mutate(

        dow.fct = factor(dow),
        week.fct = factor(week),
        year.fct = factor(year),
        zip.num2 = zip.num,

        heat.grp = inla.group(heat,100),
        heat.grp.scale = as.numeric(scale(heat.grp)),
        heat.grp.scale2 = heat.grp.scale
        
    ) %>%
    select(
        
        er,
        cz,
        n,
        doy,
        dow.fct,
        week.fct,
        year.fct,
        
        zip,
        zip.num,
        zip.num2,
        
        heat.grp,
        heat.grp.scale,
        heat.grp.scale2
                
    ) %>%
    drop_na() %>%
    arrange(zip)

temp <- dat.spatial.org %>%
    group_by(zip) %>%
    mutate(
        heat.grp.scale.mean = mean(heat.grp.scale,na.rm=TRUE)
    ) %>%
    ungroup() %>%
    mutate(
        diff.abs = abs(heat.grp.scale - heat.grp.scale.mean)
    ) %>%
    group_by(zip) %>%
    reframe(
        heat.grp.scale.compare = heat.grp.scale[which.min(diff.abs)]
    ) %>%
    ungroup() %>%
    select(
        zip,
        heat.grp.scale.compare
    )

dat.spatial <- left_join(dat.spatial.org,temp,by="zip")

preds <- dat.spatial %>%
    mutate(
        er=NA
    ) %>%
    group_by(zip) %>%
    mutate(
        heat.grp.scale = heat.grp.scale.compare, 
        heat.grp.scale2 = heat.grp.scale,
    ) %>%
    ungroup()

dat.spatial.preds <- bind_rows(dat.spatial,preds)

print("DID YOU ADD IN YEAR.FCT for f.spatial???!!!")

source("../support/priors.R")
f.inla.spatial <- er ~ 1 +
    year.fct + ## Fix this
    f(cz,
      model="iid",
      constr=TRUE) +
    f(doy,
      model="rw1",
      constr=TRUE,
      scale.model=TRUE,
      hyper = prior.list$dunif 
      ) +
    dow.fct +
    f(zip.num,
      model = "bym2",
      graph=H,
      hyper=prior.list$dunif,
      constr=TRUE,
      scale.model = TRUE,
      adjust.for.con.comp=TRUE) +
    heat.grp.scale +
    f(zip.num2,
      heat.grp.scale,
      model = "besag",
      graph=H,
      hyper=prior.list$dunif,
      constr=TRUE,
      scale.model = TRUE,
      adjust.for.con.comp=TRUE) +
    f(heat.grp.scale2,
      model="rw1",
      hyper=prior.list$dunif,
      constr=TRUE)

J <- n_distinct(dat.spatial$zip)
lcs <- inla.make.lincombs(
    heat.grp.scale = rep(1,J),
    zip.num2 = diag(J)
)
names(lcs) <- dat %>%
    select(zip) %>%
    distinct() %>%
    pull(zip)

mod <- inla( 
    f.inla.spatial,
    lincomb=lcs,
    family="poisson",
    control.predictor = list(compute=TRUE,link=1),
    control.compute=list(config=TRUE),
    ##control.compute = list(dic = TRUE, return.marginals.predictor = TRUE),
    data=dat.spatial.preds
)

betas.marg <- mod$marginals.lincomb.derived
lincomb.tibble <- tibble(
    zip = names(betas.marg),
    lincomb.slopes = unlist(lapply(betas.marg,function(x){inla.emarginal(function(k)k,x)})),
    lincomb.probs = unlist(lapply(betas.marg,function(x){1-inla.pmarginal(0,x)}))
)

N <- nrow(dat.spatial)
y.hat <- mod$summary.fitted.values[1:N,"mean"]
y.hat.preds <- mod$summary.fitted.values[N + 1:N,"mean"]

dat.spatial.post <- dat.spatial %>%
    mutate(
        y.hat,
        y.hat.preds,
        rr = y.hat / y.hat.preds,
        diff = y.hat - y.hat.preds
    ) %>%
    left_join(
        .,
        lincomb.tibble,
        by="zip"
    ) %>%
    arrange(zip,heat.grp.scale)

######################
### Threshold Maps ###
######################

mf <- function(rr.temp,heat.grp,rr.cut,slopes) {
    
    index <- rr.temp >= rr.cut 

    my.slopes <- mean(slopes) ## should all be the same per zip
    if(all(index == FALSE) | my.slopes <= 0 | is.na(my.slopes))
        return(NA)
    
    return(median(heat.grp[index]))
    ##return(min(heat.grp[index]))
}

make.thres <- function(rr.cut) {
    
    rr.thres <- dat.spatial.post %>%
        group_by(zip) %>%
        reframe(
            heat.threshold = mf(rr,heat.grp,rr.cut,lincomb.slopes) 
        ) %>%
        ungroup() 
    
    rr.thres
}
 
library(viridis)
library(RColorBrewer)

make.thres.map <- function(rr.cut.vec,file.name,min.percent=0.0,max.percent=1.0,num.cut=10) {

    ma.combo <- NULL
    for(i in 1:length(rr.cut.vec)) {
        rr.cut <- rr.cut.vec[i]
        ma.temp <- left_join(cali.shp,make.thres(rr.cut),by="zip")
        ma.temp$rr.cut <- round(rr.cut,3)
        ma.combo <- bind_rows(ma.combo,ma.temp)
    }
    
    ## min.thres <- min(ma.combo$heat.threshold,na.rm=TRUE)
    min.thres <- quantile(ma.combo$heat.threshold,prob=min.percent,na.rm=TRUE)
    ## max.thres <- max(ma.combo$heat.threshold,na.rm=TRUE)
    max.thres <- quantile(ma.combo$heat.threshold,prob=max.percent,na.rm=TRUE)
    my.breaks <- seq(min.thres,max.thres,length.out=num.cut)
    
    png(file.name, width = 1000, height = 1000)
    g <- ggplot(data = ma.combo) +
        geom_sf(aes(fill = heat.threshold),color=NA) +
        scale_fill_viridis_c(
            breaks=my.breaks,
            limits=c(min.thres,max.thres)
        ) +
        guides(fill = guide_legend(reverse = TRUE, title = "Heat Threshold")) +
        facet_wrap(~rr.cut) 
    plot(g)
    dev.off()
}

rr.cut.vec <- seq(1,1.5,length.out=4)
file.name <- "thresholds.rr.png"
make.thres.map(rr.cut.vec,file.name,min.percent=0.01,max.percent=1.0)

### check
a <- mod$summary.random$zip.num[,"mean"][dat.spatial$zip.num]
print(cor(a,dat.spatial$heat.grp))

##################
### Line Plots ###
##################

zips <- dat.spatial.post %>% distinct(zip) %>% pull(zip)
years <- dat.spatial.post %>% distinct(year.fct) %>% pull(year.fct)

for(pz in c(1,100)) {
    
    pick.zip <- pz
    pick.year <- 1
    dat.plot.org <- dat.spatial.post %>%
        filter(
            zip == zips[pick.zip],
            year.fct == years[pick.year]
        ) %>%
        arrange(doy)
    
    mm <- inla(rr ~ 1 + f(doy,model="rw2",constr=TRUE),data=dat.plot.org)
    
    dat.plot <- dat.plot.org %>%
        mutate(
            rr.mean = mm$summary.fitted.values[,"mean"],
            rr.lcl = mm$summary.fitted.values[,"0.025quant"],
            rr.ucl = mm$summary.fitted.values[,"0.975quant"],
            heat.grp.fct = factor(heat.grp,levels=unique(heat.grp))
        )
    
    file.name <- str_c("plot.",zips[pick.zip],".year.",years[pick.year],".png")
    png(file.name,width = 1000, height = 1000)
    g <- dat.plot %>%
        ggplot(aes(x = doy, y = rr)) +
        geom_point(size = 1) +
        geom_line(aes(x = doy, y = rr.mean),linetype=1,col="blue") +
        geom_line(aes(x = doy, y = rr.lcl),linetype=2,col="purple") +
        geom_line(aes(x = doy, y = rr.ucl),linetype=2,col="purple") +
        geom_hline(yintercept=rr.cut.vec,linetype=1) + #:length(rr.cut.vec)) +
        scale_x_continuous(
            breaks=dat.plot$doy[seq(1,length(dat.plot$doy),by=3)],
            labels=round(dat.plot$heat.grp[seq(1,length(dat.plot$heat.grp),by=3)],0)#,
            ##guide = guide_axis(n.dodge = 2)
        ) +
        scale_y_continuous(
            breaks=c(0,3,rr.cut.vec),
            labels=round(c(0,3,rr.cut.vec),2),
            limits=c(0.5,3)
        ) +
        theme(axis.text.x = element_text(angle = -45, hjust = 0,size=12)) +
        xlab("Heat by Day") + 
        ylab("RR") +
        theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
        )
    plot(g)
    dev.off()

    ## Heat by value ###

    file.name.by.value <- str_c("plot.by.value.",zips[pick.zip],".png")
    png(file.name.by.value, width = 1000, height = 1000)
    plot(dat.plot$heat.grp,dat.plot$rr,xlab="Heat by Value",ylab="RR",xaxt="n",ylim=c(0.5,3))
    axis(1,
         at=seq(
             min(dat.plot$heat.grp),
             max(dat.plot$heat.grp),
             1
         ),
         labels=round(
             seq(
                 min(dat.plot$heat.grp),
                 max(dat.plot$heat.grp),
                 1
             )
         ))
    abline(h=rr.cut.vec)
    dev.off()
        
    ## Zip Plot ###
    
    zip.plot <- cali.shp %>%
        mutate(
            value = if_else(zip == zips[pick.zip],1,NA)
        )

    file.name.map <- str_c("map.",zips[pick.zip],".year.",years[pick.year],".png")
    png(file.name.map, width = 1000, height = 1000)
    g <- ggplot(data = zip.plot) +
        geom_sf(aes(fill = value),color=NA) +
        scale_fill_viridis_c() +
        guides(fill="none")
    plot(g)
    dev.off()

}


