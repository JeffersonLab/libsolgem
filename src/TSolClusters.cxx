#include "TSolClusters.h"
#include <assert.h>

Double_t TSolClusters::gSplitFrac = 0.1;
UInt_t TSolClusters::gMaxClusterSize = 7; 
UInt_t TSolClusters::gMaxHits = 5;
Double_t TSolClusters::gBig = 1e9;

Int_t 
TSolClusters::MakeClusters (const Double_t stripstart, const Double_t strippitch)
{
  // Find and analyze clusters for a given plane. Clusters of active strips are considered
  // a "Hit".
  //
  // The cluster analysis is a critical part of the GEM analysis. Various
  // things can and probably need to be done right here already: splitting 
  // oversized clusters, detecting noise hits/bogus clusters, detecting and
  // fitting overlapping clusters etc. 
  //
  // This analysis may even need to be re-done after preliminary tracking to
  // see if the clustering can be improved using candidate tracks.
  // Additionally, correlated amplitude information from a second readout
  // direction in the same readout plane could be used here. These advanced
  // procedures would require significant redesign of the code: 
  // all raw strip info will have to be saved and prcessed at a later point, 
  // similar to the finding of hit pairs in like-oriented planes of the MWDC.
  //
  // For the moment, we implement a very simple algorithm: any cluster of 
  // strips larger than what a single cluster should be is assumed to be two or
  // more overlapping hits, and the cluster will be split as follows: anything
  // that looks like a local peak followed by a valley will be considered an
  // actual cluster. The parameter frac = gSplitFrac (0.0 ... 1.0) can
  // be used for some crude tuning. frac > 0.0 means that a peak is
  // only a peak if the amplitude drops below (1-frac), so
  // frac = 0.1 means: trigger on a drop below 90 % etc. Likewise for the
  // following valley: the bottom is found if the amplitude rises again
  // by (1 + frac), so frac = 0.1 means: trigger on a rise above 110 % etc.
  //
  // Strip origin and pitch are passed as arguments, so cluster position, 
  // width, and resolution have same units as the passed values.

  Double_t frac_down = 1.0 - gSplitFrac, frac_up = 1.0 + gSplitFrac;
  Cluster_t* prevHit = 0;
  UInt_t nHits = 0;

  HitMap_t::iterator next = fRawHits.begin();
  while (next != fRawHits.end()) 
    {
      HitMap_t::iterator start = next, cur = next;
      ++next;
      assert (next == fRawHits.end() or ((*next).first > (*cur).first));
      while (next != fRawHits.end() and ((*next).first - (*cur).first == 1)) 
	{
	  ++cur;
	  ++next;
	}
      
      // Now the cluster candidate is between start and cur
      assert ((*cur).first >= (*start).first);
      Int_t type = 0;
      UInt_t size = (*cur).first - (*start).first + 1;
      if (size > gMaxClusterSize) 
	{
	  Double_t maxadc = 0.0, minadc = gBig;
	  HitMap_t::iterator it = start, maxpos = start, minpos = start;
	  enum EStep 
	  {kFindMax = 1, kFindMin, kDone
	  };
	  EStep step = kFindMax;
	  while (step != kDone and it != next) 
	    {
	      Double_t adc = (*it).second;
	      switch (step) 
		{
		case kFindMax:
		  // Looking for maximum
		  if (adc > maxadc) 
		    {
		      maxpos = it;
		      maxadc = adc;
		    }
		  else if (adc < maxadc * frac_down) 
		    {
		      assert (maxadc > 0.0);
		      step = kFindMin;
		      continue;
		    }
		  break;
		case kFindMin:
		  // Looking for minimum
		  if (adc < minadc) 
		    {
		      minpos = it;
		      minadc = adc;
		    }
		  else if (adc > minadc * frac_up) 
		    {
		      assert (minadc < gBig);
		      step = kDone;
		    }
		  break;
		case kDone:
		  assert (false);
		  // should never get here
		  break;
		}
	      
	      ++it;
	    }
	  
	  if (step == kDone) 
	    {
	      // Found maximum followed by minimum
	      assert (minpos != start);
	      assert (minpos != cur);
	      assert ((*minpos).first > (*maxpos).first);
	      // Split the cluster at the position of the minimum, assuming that
	      // the strip with the minimum amplitude is shared between both clusters
	      cur = minpos;
	      next = minpos;
	      // In order not to double-count amplitude, we split the signal height
	      // of that strip evenly between the two clusters. This is a very
	      // crude way of doing what we really should be doing: "fitting" a peak 
	      // shape and using the area and centroid of the curve
	      (*minpos).second /= 2.0;
	    }
	  
	  type = step;
	  size = (*cur).first - (*start).first + 1;
	  assert ((*cur).first >= (*start).first);
	}

      // Compute weighted position average. Again, a crude (but fast) substitute
      // for fitting the centroid of the peak.
      Double_t xsum = 0.0, adcsum = 0.0;
      for (; start != next; ++start) 
	{
	  Double_t pos = stripstart + (*start).first * strippitch;
	  Double_t adc = (*start).second;
	  
	  xsum += pos * adc;
	  adcsum += adc;
	}
      
      assert (adcsum > 0.0);
      // Make a new hit
      //
      // The "type" parameter indicates the result of the custer analysis:
      // 0: clean (i.e. smaller than gMaxClusterSize, no further analysis)
      // 1: large, maximum at right edge, not split
      // 2: large, no clear minimum on the right side found, not split
      // 3: split, well-defined peak found (may still be larger than maxsize)
      //
      // The resolution (sigma) of the position measurement depends on the 
      // cluster size. In particular, if the cluster consists of only a single
      // hit, the resolution is much reduced
      Double_t resolution = 0; // not computed yet
      // resolution = fResolution;
      // if (size == 1) 
      // 	{
      // 	  resolution = 0.25*strippitch;
      // 	  // The factor of 1/2*pitch is just a guess. Since with real GEMs
      // 	  // there _should_ always be more than one strip per cluster, we must
      // 	  // assume that the other strip (s) did not fire due to inefficiency.
      // 	  // As a result, the error is bigger than it would be if only ever one
      // 	  // strip fired per hit.
      // 	  // resolution = TMath::Max (0.5*strippitch, 2.0*fResolution);
      // 	  //
      // 	} else if (size == 2) 
      // 	{
      // 	  // // Again, this is a guess, to be quantified with Monte Carlo
      // 	  // resolution = 1.2*fResolution;
      // 	}
      
      Double_t pos = xsum/adcsum;
      Cluster_t* theCluster = new Cluster_t;
      theCluster->fPos = pos;
      theCluster->fSize = size;
      theCluster->fCharge = adcsum;
      theCluster->fType = type;
      theCluster->fResolution = resolution;
      fClusters.push_back (*theCluster);
      ++nHits;

      // Ensure hits are ordered by position (should be guaranteed by std::map)
      assert (prevHit == 0 or theCluster->fPos >= prevHit->fPos);
      if (prevHit) delete prevHit;
      prevHit = theCluster;
    }
  
  
  // Negative return value indicates potential problem
  if (nHits > gMaxHits)
    nHits = -nHits;
  
  return nHits;
}
