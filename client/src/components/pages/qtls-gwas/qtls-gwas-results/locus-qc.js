import React, { useState } from 'react';
import { useDispatch, useSelector } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Zoom from '../../../controls/zoom/zoom';

export function LocusQC() {
  const {
    submitted,
    qcError,
    request,
    isLoading,
    locus_qc,
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
      {submitted && qcError && (
        <LoadingOverlay
          active={true}
          content={
            <b className="text-danger">{qcError}</b>
          }
        />
      )}
      {submitted && !qcError && !isLoading && !isLoadingQC && (
        <>
          <div style={{whiteSpace: 'pre-wrap'}}>
            {locus_qc}
          </div>
          
          <Zoom
            plotURL={`api/results/${request}/${request}_QC_overlapping.svg`}
            className="border rounded p-3 mb-2"
            maxHeight="1000px"
          />

          <Zoom
            plotURL={`api/results/${request}/${request}_QC_zscore.svg`}
            className="border rounded p-3"
            maxHeight="800px"
          />

          <Zoom
            plotURL={`api/results/${request}/${request}_QC_QTLminP.svg`}
            className="border rounded p-3"
            maxHeight="800px"
          />
        </>
      )}
    </div>
  );
}
