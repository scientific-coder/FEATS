#ifndef FEATS_SEGMENTATIONS_BOTTOM_UP_SEGMENTATOIN_HXX
#define FEATS_SEGMENTATIONS_BOTTOM_UP_SEGMENTATOIN_HXX


namespace feats{
  namespace segmentations{
    template<typename Seg, typename Nb=int> struct bottom_up_segmentation {
      typedef Seg segment_type;
      typedef Nb nb_segments_type;
    };
  }
}
#endif
