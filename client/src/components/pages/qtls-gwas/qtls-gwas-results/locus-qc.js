import React from 'react';
import { useDispatch, useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export function LocusQC() {
  const {
    submitted,
    isError,
    request,
    isLoadingQC
  } = useSelector((state) => state.qtlsGWAS);

  return (
    <div className="px-3 py-2" style={{ minHeight: '500px' }}>
      {!submitted && (
        <LoadingOverlay
          active={true}
          content={
            <>Select data in the left panel and click <b>Calculate</b> to see results here.</>
          }
        />
      )}
      {submitted && isError && (
        <LoadingOverlay
          active={true}
          content={
            <b className="text-danger">Please check input files. Reset form to try again.</b>
          }
        />
      )}
      {submitted && !isError && !isLoadingQC && (
        <>
          <img src={`api/results/${request}/${request}_QC_overlapping.svg`} style={{marginLeft: 'auto', marginRight:'auto', display:'block'}} />

          <img src={`api/results/${request}/${request}_QC_zscore.svg`} style={{marginLeft: 'auto', marginRight:'auto', display:'block'}}/>
        </>
      )}
    </div>
  );
}
