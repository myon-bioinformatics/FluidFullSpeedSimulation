import Std

/-!
Exact integer accounting behind FFSS quantity transfers.

These statements describe a mathematical model, not an extraction or proof of
the Python implementation. Python regression tests check the connection to
Network.step and DelayLine.push on concrete states. One integer unit represents
one Q48 quantum; no floating-point approximation appears in these statements.
-/

namespace FFSS

/-- The same integer removed from one stock is added to the other. -/
theorem transfer_conserves (source destination moved : Int) :
    (source - moved) + (destination + moved) = source + destination := by
  omega

/-- The transfer must be nonnegative and no larger than the available stock. -/
theorem transfer_nonnegative (source destination moved : Int)
    (destination_nonnegative : 0 ≤ destination)
    (moved_nonnegative : 0 ≤ moved) (within_stock : moved ≤ source) :
    0 ≤ source - moved ∧ 0 ≤ destination + moved := by
  omega

/-- A queue's new stock plus its output equals its old stock plus its input. -/
theorem queue_push_conserves (queued sent released : Int) :
    (queued + sent - released) + released = queued + sent := by
  omega

/--
A delayed link has two storage locations: the queue and the pending output.
The old pending output arrives at the destination; this step's queue output
becomes the new pending output. Omitting either location loses conservation.
-/
theorem delayed_transfer_conserves
    (source destination queued pending sent released : Int) :
    (source - sent) + (destination + pending) +
        (queued + sent - released) + released =
      source + destination + queued + pending := by
  omega

/-- All four resulting stocks stay nonnegative under the stated bounds. -/
theorem delayed_transfer_nonnegative
    (source destination queued pending sent released : Int)
    (destination_nonnegative : 0 ≤ destination)
    (pending_nonnegative : 0 ≤ pending)
    (released_nonnegative : 0 ≤ released)
    (within_stock : sent ≤ source)
    (within_queue : released ≤ queued + sent) :
    0 ≤ source - sent ∧ 0 ≤ destination + pending ∧
      0 ≤ queued + sent - released ∧ 0 ≤ released := by
  omega

end FFSS

-- Keep the dependencies of every checked statement visible in CI output.
#print axioms FFSS.transfer_conserves
#print axioms FFSS.transfer_nonnegative
#print axioms FFSS.queue_push_conserves
#print axioms FFSS.delayed_transfer_conserves
#print axioms FFSS.delayed_transfer_nonnegative
