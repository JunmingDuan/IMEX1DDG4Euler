#ifndef INTERVAL_VOLUME_H
#define INTERVAL_VOLUME_H
/**
 * @file interval_volume.h
 * @brief 1D cell volume
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-04
 */

double interval_volume(const std::vector<double>& v) {
  return v[1] - v[0];
}

#endif //INTERVAL_VOLUME_H

